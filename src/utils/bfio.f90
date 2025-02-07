module bfio
  use iso_fortran_env
    implicit none

    integer, parameter, public :: FILEPATHLEN = 1024

    public:: print_centered_header,print_mit_license
  contains


  
    subroutine print_centered_header(header_text, io_unit)
      !--------------------------------------------------------------------!
      !   サブルーチン: print_centered_header                              !
      !   説明: 中央揃えのヘッダーを指定の I/O ユニットに出力する        !
      !--------------------------------------------------------------------!
      integer,   optional, intent(in) :: io_unit          ! 出力先のI/Oユニット番号
      character(len=*), intent(in) :: header_text  ! ヘッダーテキスト
      integer :: line_width           ! 行の幅
  
      integer :: text_length, padding

      line_width = 70
  
      ! テキストの長さを計算
      text_length = len_trim(header_text)
  
      ! 行幅よりテキストが長い場合は、エラーメッセージを出力
      if (text_length > line_width) then
        write(io_unit, '(A)') 'Error: Header text is longer than the specified line width!'
        return
      end if
  
      ! 左側の余白を計算 (中央揃え用)
      padding = (line_width - text_length) / 2
  
      ! 上下の区切り線を出力
      if (present(io_unit)) then
        write(io_unit,'(A)') repeat('-', line_width)
        write(io_unit,'(A)') repeat(' ', padding) // trim(header_text)
        write(io_unit,'(A)') repeat('-', line_width)
      else
        write(*,'(A)') repeat('-', line_width)
        write(*,'(A)') repeat(' ', padding) // trim(header_text)
        write(*,'(A)') repeat('-', line_width)
      endif
  
    end subroutine print_centered_header




    subroutine print_mit_license(io_unit)
        !--------------------------------------------------------------------!
        !   サブルーチン: print_mit_license                                  !
        !   説明: MITライセンスを指定のI/Oユニットに表示する               !
        !--------------------------------------------------------------------!
        integer,   optional, intent(in) :: io_unit  ! 出力先のI/Oユニット番号
    
        ! ライセンスの内容を出力
        if (present(io_unit)) then
            write(io_unit, '(A)') '====================================================================='
            write(io_unit, '(A)') '                                MIT License                         '
            write(io_unit, '(A)') '====================================================================='
            write(io_unit, '(A)') 'Copyright (c) 2024 Yuki Nagai'
            write(io_unit, '(A)') ''
            write(io_unit, '(A)') 'Permission is hereby granted, free of charge, to any person'
            write(io_unit, '(A)') 'obtaining a copy of this software and associated documentation'
            write(io_unit, '(A)') 'files (the "Software"), to deal in the Software without'
            write(io_unit, '(A)') 'restriction, including without limitation the rights to use,'
            write(io_unit, '(A)') 'copy, modify, merge, publish, distribute, sublicense, and/or sell'
            write(io_unit, '(A)') 'copies of the Software, and to permit persons to whom the'
            write(io_unit, '(A)') 'Software is furnished to do so, subject to the following'
            write(io_unit, '(A)') 'conditions:'
            write(io_unit, '(A)') ''
            write(io_unit, '(A)') 'The above copyright notice and this permission notice shall be'
            write(io_unit, '(A)') 'included in all copies or substantial portions of the Software.'
            write(io_unit, '(A)') ''
            write(io_unit, '(A)') 'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,'
            write(io_unit, '(A)') 'EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES'
            write(io_unit, '(A)') 'OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND'
            write(io_unit, '(A)') 'NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT'
            write(io_unit, '(A)') 'HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,'
            write(io_unit, '(A)') 'WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING'
            write(io_unit, '(A)') 'FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR'
            write(io_unit, '(A)') 'OTHER DEALINGS IN THE SOFTWARE.'
            write(io_unit, '(A)') '====================================================================='
        else
            write(*, '(A)') '====================================================================='
            write(*, '(A)') '                                MIT License                         '
            write(*, '(A)') '====================================================================='
            write(*, '(A)') 'Copyright (c) 2024 Yuki Nagai'
            write(*, '(A)') ''
            write(*, '(A)') 'Permission is hereby granted, free of charge, to any person'
            write(*, '(A)') 'obtaining a copy of this software and associated documentation'
            write(*, '(A)') 'files (the "Software"), to deal in the Software without'
            write(*, '(A)') 'restriction, including without limitation the rights to use,'
            write(*, '(A)') 'copy, modify, merge, publish, distribute, sublicense, and/or sell'
            write(*, '(A)') 'copies of the Software, and to permit persons to whom the'
            write(*, '(A)') 'Software is furnished to do so, subject to the following'
            write(*, '(A)') 'conditions:'
            write(*, '(A)') ''
            write(*, '(A)') 'The above copyright notice and this permission notice shall be'
            write(*, '(A)') 'included in all copies or substantial portions of the Software.'
            write(*, '(A)') ''
            write(*, '(A)') 'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,'
            write(*, '(A)') 'EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES'
            write(*, '(A)') 'OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND'
            write(*, '(A)') 'NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT'
            write(*, '(A)') 'HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,'
            write(*, '(A)') 'WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING'
            write(*, '(A)') 'FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR'
            write(*, '(A)') 'OTHER DEALINGS IN THE SOFTWARE.'
            write(*, '(A)') '====================================================================='

        endif
      end subroutine print_mit_license
  
  end module bfio