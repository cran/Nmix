  module f2cio
  
  use iso_c_binding

  implicit none
  private
  public :: fopen, fclose, fread, fwrite, cseek
  
  interface

     function fopen(filename, mode) bind(C,name='fopen')
       !filename: file name to associate the file stream to 
       !mode: null-terminated character string determining file access mode 
       import
       implicit none
       type(c_ptr) fopen
       character(kind=c_char), intent(in) :: filename(*)
       character(kind=c_char), intent(in) :: mode(*)
     end function fopen

     function fclose(fp) bind(C,name='fclose')
       !fp: the file stream to close 
       import
       implicit none
       integer(c_int) fclose
       type(c_ptr), value :: fp
     end function fclose
     
     function fread(buffer,size,nbytes,fp) bind(C,name='fread')
       ! buffer: pointer to the array where the read objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be read 
       ! fp: the stream to read 
       import
       implicit none
       integer(c_int) fread
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       type(c_ptr), value :: buffer 
       type(c_ptr), value :: fp
     end function fread
     
     function cseek(fp,offset,origin) bind(C,name='fseek')
       !fp: file stream to modify 
       !offset: number of characters to shift the position relative to origin 
       !origin: position to which offset is added (SEEK_SET, SEEK_CUR, SEEK_END) 
       import
       implicit none
       integer(c_int) cseek
       type(c_ptr), value :: fp
       integer(kind=c_int64_t), value :: offset
       integer(kind=c_int), value :: origin
     end function cseek
     
     function fwrite(buffer,size,nbytes,fp) bind(C,name='fwrite')
       ! buffer: pointer to the array where the write objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be written 
       ! fp: the stream to write 
       import
       implicit none
       integer(c_int) fwrite
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       type(c_ptr), value :: buffer 
       type(c_ptr), value :: fp
     end function fwrite

  end interface
     
  end module f2cio
