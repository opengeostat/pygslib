!test program to read datamine double presition

subroutine dm2csv_ep(filein, fileout, nformat)


    ! input variables
    character(len=500), intent(in)  :: filein
    character(len=500), intent(in)  :: fileout
    character(len=500), intent(in)  :: nformat

    
    ! tmp variables
    character*4 text4
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
    
    open (unit=1, file=filein, form='unformatted', access='stream')
    open (unit=2, file=fileout)
    
    ! read the File name: 1-4 and 9-12
    read (unit=1, pos=1) text4
    write (*,'(A, $)') text4
    read (unit=1, pos=9) text4
    write (*,'(A)') text4

    ! read the Database name:     17-20 and 25-28
    read (unit=1, pos=17) text4
    write (*,'(A, $)') text4
    read (unit=1, pos=25) text4
    write (*,'(A)') text4

    
    ! read the File description : 33-36,41-44,...,185-188
    p=33
    l=1
    do i=1,20
        read (unit=1, pos=p) text4
        desc(l:l+4)= text4
        p=p+8
        l=l+4
    end do
    
    write (*,*) 'one line description: ', desc
    
    ! Numeric date coded as 10000*year + 100*month + day: 193-200
    read (unit=1, pos=193) numcode
    write (*,'(A,I10)') 'Numeric date coded: ', int(numcode)

    ! Total number of fields in the file (alpha fields are counted as the number of 4-byte blocks they occupy): 201-208 193-200
    read (unit=1, pos=201) nfield
    write (*,'(A,I10)') 'n fields: ', int(nfield)
    
    ! Number of last page in the file: 209-216
    read (unit=1, pos=209) npages
    write (*,'(A,I10)') 'Number of last page in the file: ', int(npages)
    
    ! Number of last logical data record within the last page:     217-224
    read (unit=1, pos=217) nreclastpage
    write (*,'(A,I10)') 'Logical records in last page: ', int(nreclastpage)
    
    write (*,*) ''
    write (*,*) '*****************************************'
    write (*,*) '      DD parameters                      '
    write (*,*) '*****************************************'
    
    p=225
    maxlen = 0
    do i=1,int(nfield)
        ! Field name
        read (unit=1, pos=p) text4
		varnames(i) = text4
        read (unit=1, pos=p+8) text4
        varnames(i) = varnames(i) // text4 ! names are stored in two variables
        
		write (2,'(A,A,$)', Advance = 'no') trim(varnames(i)),','
        write (*,*) 'Name : ', trim(varnames(i))
        
        ! Field type ('A ' or 'N ')
        read (unit=1, pos=p+16) text4
        vartype(i) = text4(1:2)
        write (*,*) 'Type : ', vartype(i)
        
        ! SW Stored word number
        read (unit=1, pos=p+24) li
        strewrdnum(i) = int(li)
        write (*,*) 'SW  : ', strewrdnum(i)
        
        ! Word number within field (always 1 for numeric fields, 1,2,3... for text fields)
        read (unit=1, pos=p+32) li
        numofwords(i) = int(li)
        write (*,*) 'WN  : ', numofwords(i)
        maxlen = maxlen + 1!numofwords(i)
        
        ! Default value or file constant value
        read (unit=1, pos=p+32) defaultval(i)
        write (*,*) 'DV  : ', defaultval(i)
        
        p=p+56
    end do
    
    write (2,'(A)') 'ID'
    
    
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
                if (vartype(k) == 'N') then
                    read (unit=1, pos=p) li
                    write (2,numfor) li, ','
                    p=p+8
                end if 
                ! if is string caracter
                if (vartype(k)=='A') then
					read (unit=1, pos=p) text4
					write (2,'(A,A,$)',advance='no') text4, ','
					p=p+8
                end if 
            end do
            write (2,*) line
        end do
        p= i*4096+1 !next page position
    end do
    
    !read last page
    do j=1, int(nreclastpage)
        line = line+1
        do k=1, int(nfield)
            ! if is a numeric field
            if (vartype(k) == 'N') then
                read (unit=1, pos=p) li
                write (2,numfor) li, ','
                p=p+8
            end if 
            ! if is string caracter
            if (vartype(k)=='A') then
                do l=1, numofwords(i)
                    read (unit=1, pos=p) text4
                    write (2,'(A,A,$)',advance='no') text4, ','
                    p=p+8
                end do 
            end if 
        end do
        write (2,*) line
    end do
    
    
    
    close (unit =1)
    close (unit =2)


end subroutine dm2csv_ep
