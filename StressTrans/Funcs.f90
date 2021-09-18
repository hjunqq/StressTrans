    module funcs
    contains
    function lowcase(s) result(t)
    ! return string 's' in lowercase
    character(*), intent(in)	::s
    character(len(s))	::t
    integer	::i,diff
    t = s; diff = 65-97
    do i=1,len(t)
        if(ichar(t(i:i))>65.and.ichar(t(i:i))<=90)then
            ! if uppercase, make lowercase
            t(i:i) = char(ichar(t(i:i))-diff)
        endif
    enddo
    end function
    end module