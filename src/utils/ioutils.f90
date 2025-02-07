module ioutils
    use iso_fortran_env
    implicit none

    public::file_exists,read_next_valid_line,&
            to_lowercase,extract_key_value,&
            find_keyword_in_file
    


! Interface for key-value extraction
    interface extract_key_value
        module procedure extract_key_value_int
        module procedure extract_key_value_real
        module procedure extract_key_value_string
    end interface 

    contains

    subroutine file_exists(filename, exists)
        !--------------------------------------------------------------------!
        !   サブルーチン: file_exists                                        !
        !   説明: ファイルが存在するかどうかを確認し、結果を返す            !
        !--------------------------------------------------------------------!
        character(len=*), intent(in) :: filename   ! チェックするファイル名
        logical, intent(out) :: exists            ! ファイルの存在状態
  
        ! ファイルをOPEN文でチェック
        inquire(file=filename, exist=exists)
  
    end subroutine file_exists

    subroutine read_next_valid_line(unit, line, end_of_file)
        implicit none
        integer, intent(in) :: unit        ! ファイルユニット番号
        character(len=*), intent(out) :: line  ! 読み込んだ有効な行を格納
        logical, intent(out) :: end_of_file    ! ファイルの終端に到達したか

        integer :: ios  ! I/O ステータス
        !logical :: is_valid

        end_of_file = .false.
        do
            ! 1行を読み込む
            read(unit, '(A)', iostat=ios) line

            if (ios < 0) then
                ! ファイルの終端
                end_of_file = .true.
                exit
            else if (ios > 0) then
                ! 読み取りエラー
                print *, 'Error reading from file.'
                stop
            end if

            ! 空行または文頭が#で始まる行を無視
            if (len_trim(line) == 0) cycle  ! 空行を無視
            if (line(1:1) == '#') cycle  ! #で始まる行を無視
            if (line(1:1) == '!') cycle  ! #で始まる行を無視
            if (line(1:1) == '%') cycle  ! #で始まる行を無視

            ! 有効な行が見つかった場合
            exit
        end do

    end subroutine read_next_valid_line

    subroutine find_keyword_in_file(filename, target_keyword, value)
        implicit none
        ! 引数
        character(len=*), intent(in) :: filename       ! ファイル名
        character(len=*), intent(in) :: target_keyword ! 検索するキーワード
        character(len=*), intent(out) :: value         ! キーワードに続く値
        
    
        ! ローカル変数
        logical :: found                 ! キーワードが見つかったか
        character(len=256) :: line
        character(len=100) :: keyword
        integer :: ios, unit
        logical :: file_opened
    
        found = .false.
        value = ''
        file_opened = .false.
    
        ! ファイルを開く
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
          print *, "Error: Cannot open file ", trim(filename)
          return
        end if
        file_opened = .true.
    
        ! ファイルを1行ずつ読み取り
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
    
          ! 行を分割してキーワードと値を取得
          read(line, *, iostat=ios) keyword, value
          if (ios == 0) then
            if (to_lowercase(trim(keyword)) == to_lowercase(trim(target_keyword))) then
              found = .true.
              exit
            end if
          end if
        end do
    
        ! ファイルを閉じる
        if (file_opened) close(unit)
        if (.not. found) then
            write(*,*)  target_keyword," is not found in ",trim(adjustl(filename))
            stop
        endif
    
    end subroutine find_keyword_in_file

    function to_lowercase(input_string) result(lowercase_string)
        implicit none
        character(len=*), intent(in) :: input_string
        character(len=len(input_string)) :: lowercase_string
        integer :: i, char_code
    
        ! 初期化
        lowercase_string = input_string
    
        ! 各文字を小文字に変換
        do i = 1, len_trim(input_string)
            char_code = iachar(input_string(i:i))
            ! 大文字 (A-Z) を小文字 (a-z) に変換
            if (char_code >= iachar('A') .and. char_code <= iachar('Z')) then
                lowercase_string(i:i) = achar(char_code + 32)
            end if
        end do
      end function to_lowercase


      subroutine extract_key_value_string(input_string, key, value)
        implicit none
        character(len=*), intent(in) :: input_string  ! 入力文字列 (例: "BASIS type=Chebyshev")
        character(len=*), intent(in) :: key           ! 検索するキー (例: "type")
        character(len=*), intent(out) :: value        ! キーに対応する値を格納
        integer :: pos_start, pos_end, len_input, key_len
    
        ! 初期化
        value = ''
        len_input = len_trim(input_string)
        key_len = len_trim(key)
    
        ! "key=" の開始位置を探す
        pos_start = index(input_string, trim(key) // '=')
        if (pos_start == 0) then
            print *, 'Error: "', trim(key), '=" not found in the string.'
            return
        end if
    
        ! キーの後の値の開始位置を計算
        pos_start = pos_start + key_len + 1  ! key= の長さ分だけ進める
    
        ! 値の終端位置を探す（空白または文字列の終わりまで）
        pos_end = scan(input_string(pos_start:len_input), ' ') + pos_start - 2
        if (pos_end < pos_start) then
            pos_end = len_input  ! 空白が見つからなければ文字列の終わりまで
        end if
    
        ! 値を取り出す
        value = input_string(pos_start:pos_end)
      end subroutine extract_key_value_string   
      


  ! --- 整数値を抽出するサブルーチン ---
  subroutine extract_key_value_int(line, key, value)
      implicit none
      character(len=*), intent(in) :: line  ! 入力文字列
      character(len=*), intent(in) :: key   ! 検索するキー
      integer, intent(out) :: value         ! 整数値を格納
      integer :: len_line

      integer::i,j


            ! 初期化
      value = 0.0
      len_line = len_trim(line)

      i = index(line, trim(key))
      if (i == 0) then
        return
      end if
      i = i + len_trim(key)
      j = index(line(i:len_line), '=')
      if (j == 0) then
        write(0,*) 'Error: invalid string in extract_key_value.'
        write(0,*) '       ', line
        stop
      end if
      i = i + j
      read(line(i:len_line), *) value
  end subroutine extract_key_value_int

  ! --- 実数値を抽出するサブルーチン ---
  subroutine extract_key_value_real(line, key, value)
      implicit none
      character(len=*), intent(in) :: line  ! 入力文字列
      character(len=*), intent(in) :: key   ! 検索するキー
      real(real64), intent(out) :: value            ! 実数値を格納
      integer :: len_line


      integer::i,j
      
      ! 初期化
      value = 0.0
      len_line = len_trim(line)

      i = index(line, trim(key))
      if (i == 0) then
        return
      end if
      i = i + len_trim(key)
      j = index(line(i:len_line), '=')
      if (j == 0) then
        write(0,*) 'Error: invalid string in extract_key_value.'
        write(0,*) '       ', line
        stop
      end if
      i = i + j
      read(line(i:len_line), *) value

  end subroutine extract_key_value_real      
end module ioutils