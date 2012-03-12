module Bio
  class DB
    module SAM
      module Tools
        module Bcf
          extend FFI::Library
          ffi_lib Bio::DB::SAM::Library.filename
          
          # FUNCTIONS  
          attach_function :vcf_close, [:pointer], :int                                    # BcfT*                       : 0,-1
          attach_function :vcf_open, [:string,:string], :pointer                          # filename, mode[r,w,b]       : BcfT*
          attach_function :vcf_hdr_read, [:pointer], :pointer                             # BcfT*                       : BcfHdrT*
          attach_function :vcf_read, [:pointer,:pointer,:pointer], :int                   # BcfT*, BcfHdrT*, Bcf1T      : -1,n (file offset?)
          attach_function :bcf_build_refhash, [:pointer], :pointer                        # BcfHdrT*                    : void(khash_t)*
          attach_function :bcf_str2id, [:pointer,:string], :int                           # void(hash_t)*, ref_name     : tid
          attach_function :bcf_idx_load, [:string], :pointer                              # filename                    : BcfIdxT*
          attach_function :bcf_idx_query, [:pointer,:int,:int], :uint64                   # BcfIdxT, tid, begin_pos     : offset
          attach_function :bgzf_seek, [:pointer,:uint64,:int], :int                       # BcfT->fp, offset, SEEK_SET  : 0,-1
          attach_function :bcf_idx_destroy, [:pointer], :void                             # BcfIdxT*                    : nil
          attach_function :bcf_idx_build, [:string], :int                                 # filename                    : 0,-1
          attach_function :bcf_is_indel, [:pointer], :int                                 # Bcf1T*                      : int
          attach_function :bcf_subsam, [:int, :pointer, :pointer], :int                   # n_smpl, list_id*, Bcf1T*    : 0
          attach_function :bcf_hdr_subsam, [:pointer, :int, :pointer,:pointer], :pointer  # BcfHdrT*,n_sub,char**,int*  : BcfHdrT*
          
          # STRUCTS
          class Bcf1T < FFI::Struct
            ## From bcf.h
            ## a member is said to be "derived" if its content can be derived from other members
            ## - derived info: ref, alt, flt, info, fmt (<-str), n_gi (<-fmt), n_alleles (<-alt), n_smpl (<-bcf_hdr_t::n_smpl)
            layout(
              :tid, :int32,     # refID
              :pos, :int32,     # 0-based position
              :l_str, :int32,   # length of str
              :m_str, :int32,   # allocated size of ->str
              :qual, :float,    # SNP quality
              :str, :string,    # concatenated string of variable length strings in VCF (from col.2 to col.7)
              :ref, :string,    # ''
              :alt, :string,    # ''
              :flt, :string,    # all point to ->str; no memory allocation
              :info, :string,   # ''
              :fmt, :string,    # ''
              :n_gi, :int,      # number of geno fields
              :m_gi, :int,      # allocated size of geno fields
              :gi, :pointer,    # array of geno fields
              :n_alleles, :int, # number of alleles
              :n_smpl, :int     # number of samples
            )
          end
          
          class BcfT < FFI::Struct
            layout(
            :is_vcf, :int,    # vcf vs. bcf flag
            :v, :pointer,     # auxillary data structure
            :fp, :pointer     # file handle
            )
          end
          
          class BcfIdxT < FFI::Struct
            layout(
            :n, :int32,
            :index2, :pointer
            )
          end
          
          class BcfHdrT < FFI::Struct
            layout(
            :n_ref, :int32,   # number of reference sequences
            :n_smpl, :int32,  # number of samples
            :l_nm, :int32,    # length of concatenated sequence names
            :l_smpl, :int32,  # length of concatenated sample names
            :l_txt, :int32,   # length of header text 
            :name, :string,   # concatenated sequence names
            :sname, :string,  # concatenated sample names
            :txt, :string,    # header text
            :ns, :pointer,    # array of sequence names
            :sns, :pointer)   # array of sample names
          end
          
          class BcfGinfoT < FFI::Struct
            attr_accessor :format, :length, :data
            layout(
            :fmt, :int32, # packed binary string 8bit offsets, 32bit container
            :len, :int,
            :data, :pointer
            )
            def initialize(*args)
              super(*args)
              self.unpack_fmt
              self.length = self[:len]
            end
            
            # Unpacks binary data. Sets data. See bcf.c::bcf_fmt_core
            # n - number of samples
            # x - number of alleles
            def unpack_data(n,a)
              data = self[:data]
              x = a*(a+1)/2
              arr = []
              n.times do |j|
                case format
                when 'PL'
                  d = data.get_array_of_uint8(j*x,x)
                  s = []
                  (0..(x-1)).each do |k|
                    s << d[k]
                  end
                  arr << s
                when 'DP','DV'
                  d = data.get_uint16(j*2)
                  arr << d
                when 'GQ'
                  d = data.get_uint8(j)
                  arr << d
                when 'SP'
                  d = data.get_uint32(j*4)
                  arr << d
                when 'GT'
                  y = data.get_uint8(j)
                  if((y>>7)&1 != 0)
                    arr << ['.','/','.']
                  else
                    s = []
                    s << ('0'.ord+(y>>3&7)).chr
                    s << "/|"[y>>6&1]
                    s << ('0'.ord+(y&7)).chr
                    arr << s
                  end
                when 'GL'
                  d = data.get_array_of_float(j*4*x,x)
                  s = []
                  (0..(x-1)).each do |k|
                    s << d[k]
                  end
                  arr << s
                else
                  # generic tags not supported by bcf
                  arr << "."
                end
              end
              self.data = arr
            end
            
            # Unpacks fmt binary data. Sets format. See bcf.h::bcf_str2int
            def unpack_fmt
              b = self[:fmt]
              i = b & 0x0ff
              str=''
              while(i!=0)
                str = i.chr + str
                b = b >> 8
                i = b & 0x0ff
              end
              self.format = str
            end
            
            def to_s
              "#<BcfGinfoT:#{self.object_id} @format=#{self.format} @length=#{self.length} @data=#{self.data}>"
            end
          end

        end
      end
    end
  end
end