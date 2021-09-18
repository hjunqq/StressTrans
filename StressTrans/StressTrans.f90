    program StressTrans
    use datatype
    use funcs

    implicit none

    narg = iargc()
    if (narg >0) then


        do iarg = 1,narg
            call GETARG(iarg,arg)
            if(trim(arg)=="-fpath".or.trim(arg)=="-f")then
                call getarg(iarg+1,arg)
                fpath = arg
            elseif(trim(arg)=="-point".or.trim(arg)=="-p")then
                call getarg(iarg+1,arg)
                read(arg,*)pcenter
            elseif(trim(arg)=="-angular".or.trim(arg)=="-a")then
                call getarg(iarg+1,arg)
                read(arg,*)angular
            elseif(trim(arg)=="-type".or.trim(arg)=="-t")then
                call getarg(iarg+1,arg)
                read(arg,*)resType
            elseif(trim(arg)=="-help".or.trim(arg)=="-h")then
                write(*,"(A70)")"Stress Transform Program                                             "
                write(*,"(A70)")"Copyright (C) 2021 HSTAR. All rights reserved.                       "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"Usage: StressTrans <-option> [-option-parameters]                    "
                write(*,"(A70)")"    or StressTrans                                                   "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"Available Options:                                                   "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"    -fpath, or -f      the data file path                            "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"    -point, or -p      the center point of pipe                      "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"    -angular, or -a      k                                           "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"    -type, or -t output type:                                        "
                write(*,"(A70)")"                             plain - for origin                      "
                write(*,"(A70)")"                             ascii - for gidpost                     "
                write(*,"(A70)")"                             binary - for gidpost binary             "
                write(*,"(A70)")"                             hdf5 - for gidpost hdf5                 "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"Example:                                                             "
                write(*,"(A70)")"                                                                     "
                write(*,"(A70)")"    ./StressTrans.exe -f 2duanjingli\1 -p 0.0,0.0,0.0 -a -0.001      "
                write(*,"(A70)")"                                                                     "
                stop
            endif
        enddo

    else
        inpunit = 10
        open(inpunit,file='inp')
        read(inpunit,*)text
        read(inpunit,*)fpath
        read(inpunit,*)text
        read(inpunit,*)pcenter
        read(inpunit,*)text
        read(inpunit,*)angular
        read(inpunit,*)text
        read(inpunit,*)resType
    endif

    mshunit=1
    resunit=2
    open(mshunit,FILE=fpath(1:len_trim(fpath))//".flavia.msh")
    open(resunit,FILE=fpath(1:len_trim(fpath))//".flavia.res")

    call readmsh
    call readres

    close(mshunit)
    close(resunit)

    call gtransmtrx
    call trans_stress

    select case(trim(restype))
    case("plain")
        open(resunit,FILE=fpath(1:len_trim(fpath))//".trans.flavia.res")
        call writeres
        close(resunit)
        case default
        call fwriteresbin
    end select

    contains


    !>读取网格信息
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine readmsh
    character(150)   ::orgline,text
    integer         ::idxi,idxj,igroup,icoor,ielem,idim

    isMeshGroup = .false.
    isOldFormat = .True.
    ngroup=0
    ncoor=0
    nelem=0
    ngroup=0
    !将坐标和单元的头尾节点指针置空
    NULLIFY(GroupHead)
    NULLIFY(GroupLast)
    nullify(CoorHead)
    nullify(CoorLast)
    nullify(ElemHead)
    nullify(ElemLast)
    do
        !读取一行数据
        read(mshunit,'(A150)',END=100)orgline
        orgline = adjustl(orgline)                            !for old format
        read(orgline(1:index(orgline,' ')),'(A20)')text
        !取数据的第一个单词
        select case(lowcase(trim(text)))
        case('#')
            print *,"Read A Comment"
            isOldFormat = .False.
        case('group')
            print *,"This Result is in Group"
            isMeshGroup = .true.
            isOldFormat = .False.
        case('end')
            print *,"End Group"
            goto 100
        case('mesh')
            !开始读取单元组节点信息
            ngroup=ngroup+1
            !为单元组分配内存
            ALLOCATE(PGroup)
            PGroup.index=ngroup

            if(.not.isOldFormat)then
                idxi=index(orgline,'MESH')+len("MESH")+2          !去掉双引号，所以要多加一个字符
                idxj=index(orgline,'dimension')-1-2                !去掉双引号，需要多减一个字符
                read(orgline(idxi:idxj),"(A70)")PGroup.GroupName

                idxi=index(orgline,'dimension')+len("dimension")
                idxj=index(orgline,'ElemType')-1
                read(orgline(idxi:idxj),"(I3)")PGroup.Dim

                idxi=index(orgline,'ElemType')+len("ElemType")
                idxj=index(orgline,'Nnode')-1
                read(orgline(idxi:idxj),"(A70)")PGroup.ElemType

                idxi=index(orgline,'Nnode')+len("Nnode")
                idxj=len_trim(orgline)
                read(orgline(idxi:idxj),"(I3)")PGroup.Nnode
            else
                write(PGroup.GroupName,*)ngroup

                idxi=index(orgline,'dimension')+len("dimension")
                idxj=index(orgline,'elemtype')-1
                read(orgline(idxi:idxj),*)PGroup.Dim

                idxi=index(orgline,'elemtype')+len("elemtype")
                idxj=index(orgline,'nnode')-1
                read(orgline(idxi:idxj),*)PGroup.ElemType

                idxi=index(orgline,'nnode')+len("nnode")
                idxj=len_trim(orgline)
                read(orgline(idxi:idxj),*)PGroup.Nnode

            endif




            write(*,*)"Mesh          ",PGroup.GroupName(1:len_trim(PGroup.GroupName))
            write(*,*)"Dimension     ",PGroup.Dim
            write(*,*)"ElemType      ",PGroup.ElemType(1:len_trim(PGroup.ElemType))
            write(*,*)"Nnode         ",PGroup.Nnode
            write(*,*)
            !将读取到的单元组放到链表的结尾
            PGroup.next=>NULL()
            if(ASSOCIATED(GroupLast))then
                GroupLast.next=>PGroup
                GroupLast=>PGroup
            else
                GroupHead=>PGroup
                GroupLast=>PGroup
            endif
        case("coordinates")
            print *,"Begine To Read Coordinates"
            do
                read(mshunit,'(A150)',END=100)orgline
                if(isOldFormat)orgline = adjustl(orgline)
                read(orgline(1:index(orgline,' ')),'(A20)')text
                select case(lowcase(trim(text)))
                case("end")
                    exit
                    case default
                    ncoor=ncoor+1
                    allocate(PCoor)
                    allocate(Pcoor.val(PGroup.Dim))
                    read(orgline,*)Pcoor.index,PCoor.val
                    if(ncoor/=PCoor.Index)stop 'Coor Error!'
                    PCoor.next=>Null()
                    if(associated(CoorLast))then
                        CoorLast.next=>PCoor
                        CoorLast=>PCoor
                    else
                        CoorHead=>PCoor
                        CoorLast=>PCoor
                    endif
                end select
                PGroup.Ncoor=ncoor
            enddo
        case("elements")
            print *,"Begine To Read Elements"
            ielem=0
            do
                read(mshunit,'(A150)',END=100)orgline
                if(isOldFormat)orgline = adjustl(orgline)
                read(orgline(1:index(orgline,' ')),'(A20)')text
                select case(lowcase(trim(text)))
                case("end")
                    exit
                    case default
                    nelem=nelem+1
                    ielem=ielem+1
                    allocate(PElem)
                    allocate(PElem.node(PGroup.Nnode))
                    read(orgline,*)PElem.index,PElem.node,PElem.group
                    !if(nelem/=PElem.Index)stop 'Elem Error!'
                    PElem.next=>Null()
                    if(associated(ElemLast))then
                        ElemLast.next=>PElem
                        ElemLast=>PElem
                    else
                        ElemHead=>PElem
                        ElemLast=>PElem
                    endif
                end select
                PGroup.Nelem=ielem
            enddo
            !为单元组节点中的单元列表分配内存，然后将单元列表指针指向单元信息地址
            allocate(PGroup.Elem(ielem))
            PElem=>ElemHead
            ielem=0
            do while(associated(PElem))
                if(PElem.group .eq. PGroup.index)then
                    ielem=ielem+1
                    PGroup.Elem(ielem)=PElem
                endif
                PElem=>PElem.next
            enddo
            print *,"Thera're ",ielem,"elem in this group"
        end select
    enddo

    !将链表转换成动态数组格式，以方便以后读取
    !为单元组分配内存
100 allocate(Group(ngroup))
    !从链表的头开始依次将值赋予数组
    PGroup=>GroupHead
    igroup=0
    do while(associated(PGroup))
        igroup=igroup+1
        Group(igroup)=PGroup
        PGroup=>PGroup.next
    enddo
    allocate(Coor(ncoor))
    Pcoor=>CoorHead
    icoor=0
    do while(associated(PCoor))
        icoor=icoor+1
        Coor(icoor)=PCoor
        PCoor=>PCoor.Next
    enddo
    allocate(Elem(nelem))
    PElem=>ElemHead
    ielem=0
    do while(associated(PElem))
        ielem=ielem+1
        Elem(ielem)=PElem
        PElem=>PElem.Next
    enddo
    ndim = 0
    do igroup = 1,ngroup
        ndim = max(ndim,Group(igroup).Dim)
    enddo
    end subroutine

    !>读取结果信息
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine readres
    character(300)   ::orgline,text
    integer          ::idxi,idxj,ires,icomp,icoor,istep
    real            ::curtime
    real,allocatable::time(:)

    nres = 0
    nullify(ResHead)
    nullify(ResLast)
    do
        read(resunit,'(A300)',END=100)orgline
        read(orgline(1:index(orgline,' ')),'(A20)')text
        select case(lowcase(trim(text)))
        case("ongroup")
            print *, "Group Result"
        case("end")
            goto 100
        case("result")
            nres=nres+1
            Allocate(PRes)
            PRes.index=nres
            PRes.next=>Null()
            idxi=index(orgline,'"')+1
            idxj=len_trim(orgline)
            read(orgline(idxi:idxj),'(A150)')orgline
            idxi=index(orgline,'"')-1
            read(orgline(1:idxi),'(A70)')PRes.ResName

            idxi=index(orgline,'"')+3
            idxj=len_trim(orgline)
            read(orgline(idxi:idxj),'(A150)')orgline
            idxi=index(orgline,'"')-1
            read(orgline(1:idxi),'(A70)')PRes.AnaName

            idxi=index(orgline,'"')+2
            idxj=len_trim(orgline)
            read(orgline(idxi:idxj),'(A150)')orgline
            read(orgline,*)PRes.TimeAna

            idxi=index(orgline,' ')
            idxj=len_trim(orgline)
            read(orgline(idxi:idxj),'(A150)')orgline
            read(orgline,*)PRes.ResType

            if (Pres.ResType .eq. 'Scalar')then
                PRes.nval = 1
            elseif (PRes.ResType .eq. 'Vector')then
                PRes.nval = 4
            elseif(PRes.ResType .eq. 'Matrix')then
                PRes.nval = 6
            endif
            if(associated(ResLast))then
                ResLast.next=>PRes
                ResLast=>PRes
            else
                ResHead=>PRes
                ResLast=>PRes
            endif
            write(*,*)"ResultName           ",PRes.ResName(1:len_trim(PRes.ResName))
            write(*,*)"AnalysisName         ",PRes.AnaName(1:len_trim(PRes.AnaName))
            write(*,*)"TimeAnalysis         ",PRes.TimeAna
            write(*,*)"ResultType           ",PRes.ResType(1:len_trim(PRes.ResType))
            write(*,*)
        case("componentnames")
            allocate(PRes.CompName(PRes.nval))
            idxi=index(orgline,'"')+1
            idxj=len_trim(orgline)
            read(orgline(idxi:idxj),'(A150)')orgline
            idxi=index(orgline,'"')-1
            read(orgline(1:idxi),'(A70)')PRes.CompName(1)
            do icomp=2, PRes.nval
                idxi=index(orgline,'"')+1
                idxj=len_trim(orgline)
                read(orgline(idxi:idxj),'(A150)')orgline
                idxi=index(orgline,'"')+1
                idxj=len_trim(orgline)
                read(orgline(idxi:idxj),'(A150)')orgline
                idxi=index(orgline,'"')-1
                read(orgline(1:idxi),'(A70)')PRes.CompName(icomp)
            enddo
            if(len_trim(orgline)==0)then
                PRes.nval= PRes.nval-1
            endif
        case("values")
            nullify(PRes.ValHead)
            nullify(PRes.ValLast)
            icoor=0
            do
                read(resunit,'(A150)')orgline
                read(orgline(1:index(orgline,' ')),'(A20)')text
                select case(trim(text))
                case("End")
                    exit
                    case default
                    icoor=icoor+1
                    allocate(PRes.PVal)
                    allocate(PRes.PVal.dat(PRes.nval))
                    read(orgline,*)PRes.PVal.index,PRes.PVal.dat
                    if(icoor/=PRes.PVal.index)stop 'Res Error!'
                    PRes.PVal.next=>Null()
                    if(associated(PRes.ValLast))then
                        PRes.ValLast.next=>PRes.PVal
                        PRes.ValLast=>PRes.PVal
                    else
                        PRes.ValHead=>PRes.PVal
                        PRes.ValLast=>PRes.PVal
                    endif
                end select
            enddo
            allocate(PRes.Val(icoor))
            PRes.PVal=>PRes.ValHead
            icoor=0
            do while(associated(PRes.PVal))
                icoor=icoor+1
                PRes.Val(icoor)=PRes.PVal
                PRes.PVal=>PRes.PVal.next
            enddo
            case default ! for Old Format
            if(isOldFormat)then
                nres=nres+1
                Allocate(PRes)
                PRes.index=nres
                PRes.next=>Null()
                Read(orgline,*)PRes.ResName,PRes.LoadType,PRes.TimeAna,PRes.DataType,PRes.DataLoc,PRes.DescComp

                if(Pres.DataType==1) then
                    PRes.nval = 1
                elseif(Pres.DataType==2)then
                    PRes.nval = 3
                elseif(PRes.DataType==3)then
                    PRes.nval = 6
                endif

                if(PRes.DescComp>0)then
                    Allocate(Pres.CompName(Pres.nval))
                    do icomp = 1, Pres.nval
                        read(resunit,*)Pres.CompName(icomp)
                    enddo
                endif

                write(*,*)"ResultName           ",PRes.ResName(1:len_trim(PRes.ResName))
                write(*,*)"AnalysisName         ",PRes.AnaName(1:len_trim(PRes.AnaName))
                write(*,*)"TimeAnalysis         ",PRes.TimeAna
                write(*,*)"ResultType           ",PRes.DataType
                write(*,*)

                if(associated(ResLast))then
                    ResLast.next=>PRes
                    ResLast=>PRes
                else
                    ResHead=>PRes
                    ResLast=>PRes
                endif

                nullify(PRes.ValHead)
                nullify(PRes.ValLast)
                icoor=0
                do
                    read(resunit,'(A300)')orgline
                    read(orgline(1:index(orgline,' ')),'(A20)')text
                    icoor=icoor+1
                    allocate(PRes.PVal)
                    allocate(PRes.PVal.dat(PRes.nval))
                    read(orgline,*)PRes.PVal.index,PRes.PVal.dat
                    if(icoor/=PRes.PVal.index)stop 'Res Error!'
                    PRes.PVal.next=>Null()
                    if(associated(PRes.ValLast))then
                        PRes.ValLast.next=>PRes.PVal
                        PRes.ValLast=>PRes.PVal
                    else
                        PRes.ValHead=>PRes.PVal
                        PRes.ValLast=>PRes.PVal
                    endif
                    if(icoor==ncoor)exit
                enddo
                allocate(PRes.Val(icoor))
                PRes.PVal=>PRes.ValHead
                icoor=0
                do while(associated(PRes.PVal))
                    icoor=icoor+1
                    PRes.Val(icoor)=PRes.PVal
                    PRes.PVal=>PRes.PVal.next
                enddo
            endif

        end select
    enddo
100 allocate(Res(nres))
    allocate(Time(nres))
    PRes=>ResHead
    ires=0
    do while(associated(PRes))
        ires=ires+1
        Res(ires)=PRes
        Time(ires)=Pres.TimeAna
        PRes=>Pres.Next
    enddo
    nstep = 1
    if(nres>0)then
        curtime = Time(1)
        do ires= 1,nres
            if(Time(ires)/=curtime) then
                curtime = (Time(ires))
                nstep = nstep + 1
            endif
        enddo
    endif
    end subroutine

    !>生成转换矩阵
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine gtransmtrx
    integer:: icoor
    real:: delta_coor(3), x_dir(3),rot(3,3),dy
    real:: costheta,sintheta,theta


    x_dir = (/1.0,0.0,0.0/)

    do icoor = 1, ncoor
        delta_coor = Coor(icoor)%Val - pCenter
        dy = delta_coor(2)
        delta_coor(2) = 0.0
        delta_coor(3) = delta_coor(3) + dy*angular
        theta = ACOS(DOT_PRODUCT(delta_coor,x_dir)/DOT_PRODUCT(delta_coor,delta_coor))
        theta = theta - pi/2.0
        costheta = COS(theta)
        sintheta = SIN(theta)
        rot(1,1) = costheta
        rot(1,2) = 0.
        rot(1,3) = -sintheta
        rot(2,1) = 0.
        rot(2,2) = 1.
        rot(2,3) = 0.
        rot(3,1) = sintheta
        rot(3,2) = 0.
        rot(3,3) = costheta
        allocate(Coor(icoor)%trans(3,3))
        Coor(icoor)%trans = rot
    enddo

    endsubroutine gtransmtrx

    !>生成转换矩阵
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine trans_stress
    integer:: istep, ires,icoor
    real:: stress_mat(3,3)

    allocate(ResTrans(nstep))

    istep = 0
    do ires = 1, nres
        if(Res(ires).ResName .eq. 'STRESS')then
            istep = istep +1
            restrans(istep) = res(ires)

            print *, Res(ires).TimeAna
            do icoor = 1, ncoor
                stress_mat(1,1) = res(ires)%val(icoor)%dat(1)
                stress_mat(2,2) = res(ires)%val(icoor)%dat(2)
                stress_mat(3,3) = res(ires)%val(icoor)%dat(3)
                stress_mat(1,2) = res(ires)%val(icoor)%dat(4)
                stress_mat(2,1) = res(ires)%val(icoor)%dat(4)
                stress_mat(2,3) = res(ires)%val(icoor)%dat(5)
                stress_mat(3,2) = res(ires)%val(icoor)%dat(5)
                stress_mat(1,3) = res(ires)%val(icoor)%dat(6)
                stress_mat(3,1) = res(ires)%val(icoor)%dat(6)

                stress_mat = matmul(coor(icoor)%trans,matmul(stress_mat,transpose(coor(icoor)%trans)))

                restrans(istep)%val(icoor)%dat(1) = stress_mat(1,1)
                restrans(istep)%val(icoor)%dat(2) = stress_mat(2,2)
                restrans(istep)%val(icoor)%dat(3) = stress_mat(3,3)
                restrans(istep)%val(icoor)%dat(4) = stress_mat(1,2)
                restrans(istep)%val(icoor)%dat(5) = stress_mat(2,1)
                restrans(istep)%val(icoor)%dat(6) = stress_mat(3,1)

            enddo
        endif
    enddo



    end subroutine trans_stress

    !>输出
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine writeres
    integer:: istep,icomp,icoor
    type(Resinfo),pointer::pres
    do istep = 1, nstep
        pres=>restrans(istep)
        write(resunit,101)PRes.ResName(1:len_trim(PRes.ResName)),PRes.LoadType,PRes.TimeAna,PRes.DataType,PRes.DataLoc,PRes.DescComp
        do icomp = 1, Pres.nval
            write(resunit,*)Pres.CompName(icomp)
        enddo
        do icoor = 1, ncoor
            write(resunit,100)icoor,PRes.Val(icoor).dat
        enddo

    enddo
100 format(i10,10(2x,e20.8))
101 format(a15,i8,f12.6,5i8)
    endsubroutine writeres

    subroutine writeresbin
    use gidpost
    TYPE(GiD_File) :: fdm, fdr
    integer::idx,icoor,igroup,ielem,eidx,istep
    REAL(8):: Sxx,Syy,Szz,Sxy,Syz,Sxz
    TYPE(GiD_ElementType):: etype


    select case(trim(restype))
    case('binary')
        call GiD_OpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostBinary)
    case('ascii')
        call GiD_OpenPostResultFile(fpath(1:len_trim(fpath))//'.post.msh',GiD_PostAscii)
        call GiD_OpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostAscii)
    case('hdf5')
        call GiD_OpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostHDF5)
    end select

    !fdm = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.msh',GiD_PostHDF5)
    !fdr = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostHDF5)
    !fdm = fdr
    !call GiD_OpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostHDF5)

    do igroup = 1, ngroup
        PGroup=>group(igroup)
        select case(trim(pgroup.elemtype))
        case ('Point')
            etype = Gid_Point
        case ('Linear')
            etype = Gid_Linear
        case ('Triangle')
            etype = Gid_Triangle
        case ('Quadrilateral')
            etype = Gid_Quadrilateral
        case ('Tetrahedra')
            etype = Gid_Tetrahedra
        case ('Hexahedra')
            etype = Gid_Hexahedra
        case ('Prism')
            etype = Gid_Prism
        case ('Piramid')
            etype = Gid_Piramid
        case ('Sphere')
            etype = Gid_Sphere
        case ('Circle')
            etype = Gid_Circle
        end select

        !CALL GiD_fBeginMesh(fdm,PGroup.GroupName,GiD_3D,etype,PGroup.Nnode)
        !CALL GiD_BeginMeshColor(PGroup.GroupName,GiD_3D,etype,PGroup.Nnode,0.7d0,0.7d0,0.4d0)
        CALL GiD_BeginMesh(PGroup.GroupName,GiD_3D,etype,PGroup.Nnode)
        CALL GiD_MeshUnit('m');


        !CALL GiD_fBeginCoordinates(fdm)
        CALL GiD_BeginCoordinates

        if(igroup==1)then
            do icoor =1, ncoor
                !CALL GiD_fWriteCoordinates(fdm,icoor,coor(icoor)%val(1),coor(icoor)%val(2),coor(icoor)%val(3))
                CALL GiD_WriteCoordinates(icoor,coor(icoor)%val(1),coor(icoor)%val(2),coor(icoor)%val(3))
            enddo
        endif
        !CALL GiD_fEndCoordinates(fdm)
        CALL GiD_EndCoordinates

        !CALL GiD_fBeginElements(fdm)
        CALL GiD_BeginElements

        do ielem = 1,pGroup.nelem
            !CALL GiD_fWriteElement(fdm,pgroup.elem(ielem).index,pgroup.elem(ielem).node)
            CALL GiD_WriteElement(pgroup.elem(ielem).index,pgroup.elem(ielem).node)

        enddo
        !CALL GiD_fEndElements(fdm)
        CALL GiD_EndElements

        !CALL GiD_fEndMesh(fdm)
        CALL GiD_EndMesh

    enddo

    do istep=1,nstep
        pres=>restrans(istep)
        !CALL GiD_fBeginResultHeader(fdr,'Stress','Analysis',1.0d0,GiD_Matrix,GiD_onNodes,GiD_NULL)
        CALL GiD_BeginResultHeader('Stress','Analysis',pres.TimeAna,GiD_Matrix,GiD_onNodes,GiD_NULL)
        CALL GiD_ResultValues

        do icoor = 1, ncoor
            Sxx = PRes.Val(icoor).dat(1)
            Syy = PRes.Val(icoor).dat(2)
            Szz = PRes.Val(icoor).dat(3)
            Sxy = PRes.Val(icoor).dat(4)
            Syz = PRes.Val(icoor).dat(5)
            Sxz = PRes.Val(icoor).dat(6)
            !call GiD_fWrite3DMatrix(fdr,icoor,Sxx,Syy,Szz,Sxy,Syz,Sxz)

            call GiD_Write3DMatrix(icoor,Sxx,Syy,Szz,Sxy,Syz,Sxz)
        enddo
        CALL GiD_EndResult
    enddo

    CALL GiD_ClosePostResultFile

    end subroutine writeresbin
    
    subroutine fwriteresbin
    use gidpost
    TYPE(GiD_File) :: fdm, fdr
    integer::idx,icoor,igroup,ielem,eidx,istep
    REAL(8):: Sxx,Syy,Szz,Sxy,Syz,Sxz
    TYPE(GiD_ElementType):: etype
    real(8)::rgb(3)


    select case(trim(restype))
    case('binary')
        fdr = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostBinary)
        fdm = fdr
    case('ascii')
        fdm = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.msh',GiD_PostAscii)
        fdr = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostAscii)
    case('hdf5')
        fdr = GiD_fOpenPostResultFile(fpath(1:len_trim(fpath))//'.post.res',GiD_PostHDF5)
        fdm = fdr
    end select

    do igroup = 1, ngroup
        PGroup=>group(igroup)
        select case(trim(pgroup.elemtype))
        case ('Point')
            etype = Gid_Point
        case ('Linear')
            etype = Gid_Linear
        case ('Triangle')
            etype = Gid_Triangle
        case ('Quadrilateral')
            etype = Gid_Quadrilateral
        case ('Tetrahedra')
            etype = Gid_Tetrahedra
        case ('Hexahedra')
            etype = Gid_Hexahedra
        case ('Prism')
            etype = Gid_Prism
        case ('Piramid')
            etype = Gid_Piramid
        case ('Sphere')
            etype = Gid_Sphere
        case ('Circle')
            etype = Gid_Circle
        end select
        
        call RANDOM_SEED()
        call RANDOM_NUMBER(rgb)

        CALL GiD_fBeginMeshColor(fdm,PGroup.GroupName,GiD_3D,etype,PGroup.Nnode,rgb(1),rgb(2),rgb(3))



        CALL GiD_fBeginCoordinates(fdm)

        if(igroup==1)then
            do icoor =1, ncoor
                CALL GiD_fWriteCoordinates(fdm,icoor,coor(icoor)%val(1),coor(icoor)%val(2),coor(icoor)%val(3))
            enddo
        endif
        CALL GiD_fEndCoordinates(fdm)

        CALL GiD_fBeginElements(fdm)

        do ielem = 1,pGroup.nelem
            CALL GiD_fWriteElement(fdm,pgroup.elem(ielem).index,pgroup.elem(ielem).node)

        enddo
        CALL GiD_fEndElements(fdm)

        CALL GiD_fEndMesh(fdm)

    enddo

    do istep=1,nstep
        pres=>restrans(istep)
        CALL GiD_fBeginResultHeader(fdr,'Stress','Analysis',pres.TimeAna,GiD_Matrix,GiD_onNodes,GiD_NULL)
        CALL GiD_fResultValues(fdr)

        do icoor = 1, ncoor
            Sxx = PRes.Val(icoor).dat(1)
            Syy = PRes.Val(icoor).dat(2)
            Szz = PRes.Val(icoor).dat(3)
            Sxy = PRes.Val(icoor).dat(4)
            Syz = PRes.Val(icoor).dat(5)
            Sxz = PRes.Val(icoor).dat(6)
            call GiD_fWrite3DMatrix(fdr,icoor,Sxx,Syy,Szz,Sxy,Syz,Sxz)

        enddo
        CALL GiD_fEndResult(fdr)
    enddo

    CALL GiD_fClosePostResultFile(fdr)

    end subroutine fwriteresbin

    end program