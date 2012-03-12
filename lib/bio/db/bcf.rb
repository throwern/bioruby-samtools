require 'bio/db/sam/library'
require 'bio/db/sam/bcf'

# module LibC
#   extend FFI::Library
#   ffi_lib FFI::Library::LIBC
#   attach_function :free, [ :pointer ], :void
#   # call #attach_function to attach to malloc, free, memcpy, bcopy, etc.
# end

module Bio
  class DB
    module SAM
      class Bcf
        include Bio::DB::SAM::Tools::Bcf
        
        attr_reader :bcf_name, :samples, :sequence, :txt
        
        # constructor
        def initialize(filename=nil)
          # open file
          unless File.exist?(filename)
            raise "File not Found"
          end 
          # open file
          @bcf_name = File.expand_path(filename)
          open
          # parse header attrs
          init_seq
          init_samples          
          @txt = header[:txt]
        end
        
        # CLASS METHODS
        def self.vcf_to_bcf(vcf_file,seq_file_name=nil)
          #check if file exists
          Dir.chdir(File.dirname(vcf_file))
          vcf_name = File.basename vcf_file
          unless(File.exist?(vcf_name))
            puts "File not found"
            return nil
          end
          # get the file base
          file_base = vcf_name.chomp(File.extname(vcf_name))
          # parse samples if not given
          unless(seq_file_name)
            puts "No sequence file supplied, trying to parse from file..."
            vcf_in = File.open(vcf_name,"r")
            seq_file_name = file_base+".chrom"
            seq_out = File.open(seq_file_name,"w")
            sequence_count = 0
            line = ''
            while ! line.match(/^#CHROM/)
              line = vcf_in.readline
              if(m = line.match(/##contig=<ID=(.+),length=(\d+)>/))
                seq_out.puts m[1]
                sequence_count+=1
              end
            end
            seq_out.close
            vcf_in.close
            if(sequence_count > 0)
              puts "Found #{sequence_count} sequences"
            else
              puts "None found... Please supply a list of sequence names"
              return nil
            end
          end
          # run bfctools external to create the new bcf file
          puts "launching external converter this could take a while..."
          comm = "#{File.join(File.expand_path(File.dirname(__FILE__)),'sam','external','bcftools')}"
          stdin, stdout, stderr = Open3.popen3("#{comm} view -bSu -D #{seq_file_name} #{vcf_name} > #{file_base}.bcf")
          # check for errors and return
          e = stderr.collect(&:to_s).reject{|l| l.nil? || l.empty?}
          s = stdout.collect(&:to_s)
          puts s.join("\n")
          unless e.empty?
            puts e.join("\n")
            puts "error"
            return nil
          else
            puts "Done"
            return File.expand_path(file_base+".bcf")
          end
        end
        
        # INSTANCE METHODs
        
        # create file handles
        def open
          # check file
          begin
            @bcf_file_p ||= vcf_open(self.bcf_name,"rb")
          rescue => e
            puts "Error opening bcf file #{self.bcf_name}\n#{e.message}"
            return nil
          end
          # check index - create if missing
          begin
            @bcf_index_p ||= bcf_idx_load(self.bcf_name)
            if(bcf_index_p.null?)
              puts "No index found, attempting to build..."
              if(bcf_idx_build(self.bcf_name))
                puts "done"
                @bcf_index_p = bcf_idx_load(self.bcf_name)
              else
                return nil
              end
            end
          rescue => e
            puts "Error with index file \n #{e.message}"
            return nil
          end
          return self
        end
        
        # destroy memory structures
        def close
          bcf_idx_destroy(@bcf_index_p) if @bcf_index_p
          vcf_close(@bcf_file_p) if @bcf_file_p
          @bcf_file = nil
          @bcf_file_p = nil
          @bcf_index_p = nil
          @header_p = nil
          @str2id_p = nil
        end
        
        def init_seq
          @sequence = []
          ns = header[:ns].read_pointer
          off = 0
          header[:n_ref].times do
            s = ns.get_string(off)
            @sequence << s
            off += s.length + 1
          end
        end
        
        def init_samples
          @samples = []
          sns = header[:sns].read_pointer
          off = 0
          header[:n_smpl].times do
            s = sns.get_string(off)
            @samples << s
            off += s.length + 1
          end
        end
        
        # seek bcf_file to a sequence position
        def set_position(seq, pos)
          # get Seq ID
          unless((tid = bcf_str2id(str2id_p,seq)) >= 0)
            puts "No Index found for #{seq}"
            return nil
          end
          # get region offset
          if ((off = bcf_idx_query(bcf_index_p,tid,pos)) == 0)
            puts "No records for region ..."
            return nil
          end          
          # move the file pointer
          bgzf_seek(bcf_file[:fp], off, IO::SEEK_SET)          
          return tid
        end
        
        def get_data(sequence,qstart,qend,opts={})
          sample_idx = samples.find_index(opts[:sample])
          split_hets = opts[:split_hets]
          only_variants = opts[:only_variants]          
          variants = []          
          # basic fetch function - expects pointers to Bcf1T and BcfHdrT
          fetch_function = Proc.new do |bcf_p, hdr_p|
            b_struct = Bcf1T.new(bcf_p)
            if(b_struct[:n_gi] > 0)
              if(sample_idx)
                # skip sample with no data
                next if (b_struct[:gi].get_pointer(8).get_uint8(sample_idx) >> 7 & 1) > 0
              end
              if(only_variants)
                # test all gt tags. Need non-zero for alternate match
                next if b_struct[:gi].get_pointer(8).read_array_of_uint8(b_struct[:n_smpl]).collect{|g| g & 63}.uniq == [0]
              end
            end
            # store the variant(s)
            v = Variant.new(bcf_p,hdr_p)
            variants << v
            # split heterozygous
            if(split_hets && sample_idx)
              gt = v.geno_fields.find{|g| g.format=='GT'}.data[sample_idx]
              if((gt[0]!='0'&&gt[2]=='0')||(gt[0]=='0'&&gt[2]!='0'))
                v2 = Variant.new(bcf_p,hdr_p)
                v2.variant_type='Match'
                v2.alt = "."
                variants << v2
              end
            end
          end
          # run the fetch
          fetch_with_function_raw(sequence,qstart,qend,fetch_function)
          return variants
        end
        
        # fetch data from current position
        def fetch_with_function_raw(sequence,qstart,qend,fetchFunc)
          # set sequence and position
          return nil unless (tid = set_position(sequence,qstart))
          # allocate memory
          b = FFI::MemoryPointer.new(Bcf1T,1)
          b_struct = Bcf1T.new(b)
          # main loop
          while vcf_read(bcf_file_p,header_p,b) > 0
           # check boundaries
           l = b_struct[:pos]+(b_struct[:ref]||1).length
           break if b_struct[:tid] != tid || b_struct[:pos] > qend
           next unless l >= qstart && qend > b_struct[:pos]
           # call the processing function
           fetchFunc.call(b,header_p)
          end
        end
        
        protected
        
        # bcf_t data structure
        def bcf_file
          @bcf_file ||= BcfT.new(bcf_file_p)
        end
        
        # pointer to bcf_t
        def bcf_file_p
          @bcf_file_p || (
            raise "Bcf File not loaded"
          )
        end
        
        # pointer to bcf_idx_t
        def bcf_index_p
          @bcf_index_p || (
            raise "Index File not loaded"
          )
        end
        
        def header
          @hdr ||= BcfHdrT.new(header_p)
        end
        
        # pointer to bcf_hdr_t
        def header_p
          @header_p ||= vcf_hdr_read(bcf_file_p)
        end
        
        # pointer to void(khash_t)
        def str2id_p
          @str2id_p ||= bcf_build_refhash(header_p)
        end

      end
      
      class Variant
        include Bio::DB::SAM::Tools::Bcf
        attr_accessor :is_indel, :info_tags, :geno_fields, :variant_type, :samples
        # convert struct attributes to methods
        # - gi is not included, it is overridden when re-using the C memory address and would be invalid
        @@attributes = :tid,:pos,:l_str,:m_str,:qual,:str,:ref,:alt,:flt,:info,:fmt,:n_gi,:m_gi,:n_alleles,:n_smpl
        @@attributes.each{|a| attr_accessor a}
        
        # initalize with a Bcf1T pointer, and bcfHdrT pointer
        def initialize(bcf1t_p,hdr_p)
          #header struct
          hdr = BcfHdrT.new(hdr_p)
          
          # store all of the struct values
          b_struct = Bcf1T.new(bcf1t_p)
          @@attributes.each do |a|
            self.send(a.to_s+"=",b_struct[a])
          end
          # extract genotype data
          @geno_fields = []
          self.n_gi.times do |i|
            g = BcfGinfoT.new(b_struct[:gi]+i*BcfGinfoT.size)
            g.unpack_data(hdr[:n_smpl],self.n_alleles)
            self.geno_fields << g
          end
          # misc data
          self.pos +=1
          self.qual = self.qual.to_i
          self.is_indel = bcf_is_indel(bcf1t_p) > 0
          # if self.is_indel
          #   self.variant_type = "Indel"
          # else
            if(ref.length == 1)
              if alt.length == 1
                if(alt==".")
                  self.variant_type = "Match"
                else
                  self.variant_type = "Snp"
                end
              elsif alt.length > 1
                self.variant_type = "Insertion"
              else
                self.variant_type = "Deletion"
              end
            elsif(ref.length < alt.length)
              self.variant_type ="Indel"
            else
              self.variant_type ="Deletion"
            end
          # end
          self.info_tags = Hash[ self.info.split(";").collect{|f|f.split("=")} ]
        end        
      end
      
    end
  end
end
