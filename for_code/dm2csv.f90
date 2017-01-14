!test program to read datamine double presition

subroutine dm2csv_ep(filein, fileout, nformat)


    ! input variables
    character(len=500), intent(in)  :: filein
    character(len=500), intent(in)  :: fileout
    character(len=500), intent(in)  :: nformat

    
    ! tmp variables
    character*4 text4, text4b
    character*8 text8, text8b
    character*16 text16
    real*8 li
    integer p, l    
    
    ! parameters
    character*80 desc
    real*8  numcode, nfield, npages, maxlen, nreclastpage
    character*16 varnames(500)
    character*1 vartype(500)
    integer strewrdnum(500)
    integer numofwords (500)
    real*8 defaultval (500)
    integer nrecperpage
    
    character*500 numfor
    
    ! prepare print format
    numfor = '('//trim(nformat)//',A,$)'
    
    open (unit=1, file=filein, form='unformatted', access='stream', status='old')
    open (unit=2, file=fileout)
    
    ! read the File name: 1-4 and 9-12
    read (unit=1, pos=1,err=99) text4
    read (unit=1, pos=9,err=99) text4b
    write (*,*) 'File name:', text4//text4b

    
    ! read the Database name:     17-20 and 25-28
    read (unit=1, pos=17,err=99) text4
    read (unit=1, pos=25,err=99) text4b
    write (*,*) 'Database name:', text4//text4b

    
    ! read the File description : 33-36,41-44,...,185-188
    p=33
    l=1
    desc= '' 
    do i=1,20
        read (unit=1, pos=p,err=99) text4
        desc(l:l+4)= text4
        p=p+8
        l=l+4
    end do
    write (*,*) 'File description: ', desc
    
       
    ! Numeric date coded as 10000*year + 100*month + day: 193-200
    read (unit=1, pos=193,err=99) numcode
    write (*,*) 'Numeric date coded: ', int(numcode)

    ! Total number of fields in the file (alpha fields are counted as the number of 4-byte blocks they occupy): 201-208 193-200
    read (unit=1, pos=201,err=99) nfield
    write (*,'(A,I10)') 'Number of fields: ', int(nfield)
    
    ! Number of last page in the file: 209-216
    read (unit=1, pos=209,err=99) npages
    write (*,'(A,I10)') 'Number of last page in the file: ', int(npages)
    
    ! Number of last logical data record within the last page:     217-224
    read (unit=1, pos=217,err=99) nreclastpage
    write (*,'(A,I10)') 'Logical records in last page: ', int(nreclastpage)
    
    
    write (*,*) ''
    write (*,*) '*****************************************'
    write (*,*) '      DD parameters                      '
    write (*,*) '*****************************************'
    
    p=225
    maxlen = 0
    do i=1,int(nfield)
        ! Field name
        read (unit=1, pos=p,err=99) text4
        read (unit=1, pos=p+8,err=99) text4b
        varnames(i) = text4 // text4b  ! names are stored in two variables
        
        
        write (*,*) 'Name : ', trim(varnames(i))
        
        ! Field type ('A ' or 'N ')
        read (unit=1, pos=p+16,err=99) text4
        vartype(i) = text4(1:2)
        write (*,*) 'Type : ', vartype(i)
        
        ! SW Stored word number
        read (unit=1, pos=p+24,err=99) li
        strewrdnum(i) = int(li)
        write (*,*) 'SW  : ', strewrdnum(i)
        
        ! Word number within field (always 1 for numeric fields, 1,2,3... for text fields)
        read (unit=1, pos=p+32,err=99) li
        numofwords(i) = int(li)
        write (*,*) 'WN  : ', numofwords(i)
        if (strewrdnum(i)>0) maxlen = maxlen + 1
        
        ! Default value or file constant value
        read (unit=1, pos=p+48,err=99) defaultval(i)
        write (*,*) 'DV  : ', defaultval(i)
        
        
        if (vartype(i) == 'N') then
            write (2,'(A,A,$)', Advance = 'no',err=99) trim(varnames(i)),','
        else
            if (.not. (i>1 .and. varnames(i-1)==varnames(i) .and. numofwords(i-1)==numofwords(i)-1)) then
                write (2,'(A,A,$)', Advance = 'no',err=99) trim(varnames(i)),','
            end if
        end if
        
        p=p+56
    end do
    
    write (2,'(A)',err=99) 'ID'
    
    
    nrecperpage = INT(508/maxlen)
    write (*,*) 'maxlen: ', maxlen
    write (*,*) 'number of logical data records per page: ', nrecperpage

    write (*,*) '*****************************************'
    
    ! read data
    
    
    p=4097
    line = 0
    ! read each data page 
    do i=2, int(npages-1)
        !read each record on page 
        do j=1, nrecperpage
            line = line+1
            do k=1, int(nfield)
                ! if is a numeric field
                if (vartype(k) == 'N' .and. strewrdnum(k)>0) then
                    read (unit=1, pos=p,err=99) li
                    if (li>-9.9E29 .and. li<9.9E29) then 
                        write (2,numfor,err=99) li, ','
                    else 
                        if (li<=-9.9E29) write (2,'(A,A,$)',err=99) '   ' , ','
                        if (li>= 9.9E29) write (2,'(A,A,$)',err=99) 'inf' , ','
                    end if
                    p=p+8
                end if
                if (vartype(k) == 'N' .and. strewrdnum(k)==0) then
                    write (2,numfor,err=99) defaultval(k), ','
                end if 
                ! if is string caracter
                if (vartype(k)=='A') then
                    read (unit=1, pos=p,err=99) text4
                    if(k<nfield .and. varnames(k)==varnames(k+1) .and. numofwords(k)==numofwords(k+1)-1) then
                        write (2,'(A,$)',advance='no',err=99) text4
                    else
                        write (2,'(A,A,$)',advance='no',err=99) text4, ','
                    end if 
                    p=p+8
                end if 
            end do
            write (2,*,err=99) line
        end do
        p= i*4096+1 !next page position
    end do
    
    
    !read last page
    do j=1, int(nreclastpage)
        line = line+1
        do k=1, int(nfield)
            ! if is a numeric field
            if (vartype(k) == 'N' .and. strewrdnum(k)>0) then
                read (unit=1, pos=p,err=99) li
                if (li>-9.9E29 .and. li<9.9E29) then
                    write (2,numfor,err=99) li, ','
                else 
                    if (li<=-9.9E29) write (2,'(A,A,$)',err=99) '   ' , ','
                    if (li>= 9.9E29) write (2,'(A,A,$)',err=99) 'inf' , ','
                end if
                p=p+8
            end if 
            if (vartype(k) == 'N' .and. strewrdnum(k)==0) then
                write (2,numfor,err=99) defaultval(k), ','
            end if 
            ! if is string caracter
            if (vartype(k)=='A') then
                read (unit=1, pos=p,err=99) text4
                if(k<nfield .and. varnames(k)==varnames(k+1) .and. numofwords(k)==numofwords(k+1)-1) then
                    write (2,'(A,$)',advance='no',err=99) text4
                else
                    write (2,'(A,A,$)',advance='no',err=99) text4, ','
                end if 
                p=p+8
            end if  
        end do
        write (2,*) line
    end do
    
    
    
    close (unit =1)
    close (unit =2)

    return
    
    99 error=99 !'error in data file!'

    write (*,*) 'I/O error' 
    
    close (unit =1)
    close (unit =2)
    
end subroutine dm2csv_ep



subroutine dm2csv_sp(filein, fileout, nformat)


    ! input variables
    character(len=500), intent(in)  :: filein
    character(len=500), intent(in)  :: fileout
    character(len=500), intent(in)  :: nformat

    
    ! tmp variables
    character*4 text4, text4b
    character*8 text8, text8b
    character*16 text16
    real li
    integer p, l    
    
    ! parameters
    character*80 desc
    real  numcode, nfield, npages, maxlen, nreclastpage
    character*8 varnames(500)
    character*1 vartype(500)
    integer strewrdnum(500)
    integer numofwords (500)
    real defaultval (500)
    integer nrecperpage
    
    character*500 numfor
    
    ! prepare print format
    numfor = '('//trim(nformat)//',A,$)'
    
    open (unit=1, file=filein, form='unformatted', access='stream', status='old')
    open (unit=2, file=fileout)
    
    ! read the File name: 1-4 and 9-12
    read (unit=1, pos=1,err=99) text8
    write (*,*) 'File name:', text8

    ! read the Database name:     9-6
    read (unit=1, pos=9,err=99) text8
    write (*,*) 'Database name:', text8
    
    ! read the File description : 17-96
    read (unit=1, pos=17,err=99) desc
    write (*,*) 'File description: ', desc
    
    ! Numeric date coded as 10000*year + 100*month + day: 97-100
    read (unit=1, pos=97,err=99) numcode
    write (*,*) 'Numeric date coded: ', int(numcode)

    ! Total number of fields in the file (alpha fields are counted as the number of 4-byte blocks they occupy): 201-208 193-200
    read (unit=1, pos=101,err=99) nfield
    write (*,'(A,I10)') 'Number of fields: ', int(nfield)
    
    ! Number of last page in the file: 209-216
    read (unit=1, pos=105,err=99) npages
    write (*,'(A,I10)') 'Number of last page in the file: ', int(npages)
    
    ! Number of last logical data record within the last page:     217-224
    read (unit=1, pos=109,err=99) nreclastpage
    write (*,'(A,I10)') 'Logical records in last page: ', int(nreclastpage)
    
    
    write (*,*) ''
    write (*,*) '*****************************************'
    write (*,*) '      DD parameters                      '
    write (*,*) '*****************************************'
    
    p=113
    maxlen = 0
    do i=1,int(nfield)
        ! Field name
        read (unit=1, pos=p,err=99) varnames(i)        
        write (*,*) 'Name : ', trim(varnames(i))
        
        ! Field type ('A ' or 'N ')
        read (unit=1, pos=p+8,err=99) text4
        vartype(i) = text4(1:2)
        write (*,*) 'Type : ', vartype(i)
        
        ! SW Stored word number
        read (unit=1, pos=p+12,err=99) li
        strewrdnum(i) = int(li)
        write (*,*) 'SW  : ', strewrdnum(i)
        
        ! Word number within field (always 1 for numeric fields, 1,2,3... for text fields)
        read (unit=1, pos=p+16,err=99) li
        numofwords(i) = int(li)
        write (*,*) 'WN  : ', numofwords(i)
        if (strewrdnum(i)>0) maxlen = maxlen + 1
        
        ! Default value or file constant value
        read (unit=1, pos=p+24,err=99) defaultval(i)
        write (*,*) 'DV  : ', defaultval(i)
		
        if (vartype(i) == 'N') then
            write (2,'(A,A,$)', Advance = 'no',err=99) trim(varnames(i)),','
        else
            if (.not. (i>1 .and. varnames(i-1)==varnames(i) .and. numofwords(i-1)==numofwords(i)-1)) then
                write (2,'(A,A,$)', Advance = 'no',err=99) trim(varnames(i)),','
            end if
        end if
        
        p=p+28
    end do
    
    write (2,'(A)',err=99) 'ID'
    
    
    nrecperpage = INT(508/maxlen)
    write (*,*) 'maxlen: ', maxlen
    write (*,*) 'number of logical data records per page: ', nrecperpage

    write (*,*) '*****************************************'
    
    ! read data
    
    
    p=2049
    line = 0
    ! read each data page 
    do i=2, int(npages-1)
        !read each record on page 
        do j=1, nrecperpage
            line = line+1
            do k=1, int(nfield)
                ! if is a numeric field
                if (vartype(k) == 'N' .and. strewrdnum(k)>0) then
                    read (unit=1, pos=p,err=99) li
                    if (li>-9.9E29 .and. li<9.9E29) then 
                        write (2,numfor,err=99) li, ','
                    else 
                        if (li<=-9.9E29) write (2,'(A,A,$)',err=99) '   ' , ','
                        if (li>= 9.9E29) write (2,'(A,A,$)',err=99) 'inf' , ','
                    end if
                    p=p+4
                end if
                if (vartype(k) == 'N' .and. strewrdnum(k)==0) then
                    write (2,numfor,err=99) defaultval(k), ','
                end if 
                ! if is string caracter
                if (vartype(k)=='A') then
                    read (unit=1, pos=p,err=99) text4
                    if(k<nfield .and. varnames(k)==varnames(k+1) .and. numofwords(k)==numofwords(k+1)-1) then
                        write (2,'(A,$)',advance='no',err=99) text4
                    else
                        write (2,'(A,A,$)',advance='no',err=99) text4, ','
                    end if 
                    p=p+4
                end if 
            end do
            write (2,*,err=99) line
        end do
        p= i*2048+1 !next page position
    end do
    
    
    !read last page
    do j=1, int(nreclastpage)
        line = line+1
        do k=1, int(nfield)
            ! if is a numeric field
            if (vartype(k) == 'N' .and. strewrdnum(k)>0) then
                read (unit=1, pos=p,err=99) li
                if (li>-9.9E29 .and. li<9.9E29) then
                    write (2,numfor,err=99) li, ','
                else 
                    if (li<=-9.9E29) write (2,'(A,A,$)',err=99) '   ' , ','
                    if (li>= 9.9E29) write (2,'(A,A,$)',err=99) 'inf' , ','
                end if
                p=p+4
            end if 
            if (vartype(k) == 'N' .and. strewrdnum(k)==0) then
                write (2,numfor,err=99) defaultval(k), ','
            end if 
            ! if is string caracter
            if (vartype(k)=='A') then
                read (unit=1, pos=p,err=99) text4
                if(k<nfield .and. varnames(k)==varnames(k+1) .and. numofwords(k)==numofwords(k+1)-1) then
                    write (2,'(A,$)',advance='no',err=99) text4
                else
                    write (2,'(A,A,$)',advance='no',err=99) text4, ','
                end if 
                p=p+4
            end if  
        end do
        write (2,*) line
    end do
    
    
    
    close (unit =1)
    close (unit =2)

    return
    
    99 error=99 !'error in data file!'

    write (*,*) 'I/O error' 
    
    close (unit =1)
    close (unit =2)
    
end subroutine dm2csv_sp
