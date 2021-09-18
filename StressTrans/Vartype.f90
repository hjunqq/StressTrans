    !>ȫ�ֱ���ģ�飬�����˳�������Ҫ��ȫ�ֱ���
    module datatype
    !> ��������ڵ�.
    !!@param index ���������ֵ
    !!@param val �����ֵ
    !!@param next ָ����һ������ڵ�
    type coorinfo
        integer         ::index
        real(8),allocatable::val(:)
        real,allocatable::trans(:,:)
        type(coorinfo),pointer::next
    end type
    !> ��Ԫ����ڵ�
    !!@param index ��Ԫ����ֵ
    !!@param group ��Ԫ������
    !!@param node ��Ԫ�ڵ�
    !!@param next ָ����һ����Ԫ�ڵ�
    type eleminfo
        integer         ::index, group
        integer,allocatable::node(:)
        type(eleminfo),pointer::next
    end type
    !> ��Ԫ������ڵ�
    !!@param MeshName ��Ԫ�������
    !!@param ElemType ��Ԫ������
    !!@param index ��Ԫ������ֵ
    !!@param Dim ��Ԫ��ά��
    !!@param Nnode ÿ����Ԫ�ڵ���
    !!@param Ncoor ��Ԫ���еĽڵ���
    !!@param Nelem ��Ԫ���еĵ�Ԫ��
    !!@param Next ָ����һ���ڵ�
    !!@param Elem ��Ԫ���еĵ�Ԫ�б�
    type Groupinfo
        character(70)   ::GroupName,ElemType
        integer         ::index,Dim,Nnode,Ncoor,Nelem
        type(Groupinfo),pointer::next
        type(eleminfo),dimension(:),pointer::Elem
        type(coorinfo),dimension(:),pointer::Coor
    end type
    !>���ֵ����ڵ�
    !!@param index ���ֵ����
    !!@param dat ���ֵ
    !!@param next ��һ���ڵ�
    type ResValinfo
        integer         ::index
        real(8),allocatable    ::dat(:)
        type(ResValinfo),pointer::next
    end type
    !>���������ڵ�
    !!@param index ���������
    !!@param ResName ���������
    !!@param AnaName ��������
    !!@param ResType ������ͣ���������������
    !!@param CompName ����б�ÿһ�е�����
    !!@param Val ���ֵ
    !!@param ValHead ���ֵ����ͷ�ڵ�
    !!@param ValLast ���ֵ����β�ڵ�
    !!@param PVal ��ǰ���ֵ�ڵ�
    !!@param next ��һ�������
    type Resinfo
        character(70)   ::ResName,AnaName,ResType
        character(70),allocatable::CompName(:)
        character(70)   ::LoadDesc
        integer         ::index,nval,LoadType,DataType,DataLoc,DescComp,GaussPoint
        real(8)            ::TimeAna
        type(ResValinfo),dimension(:),pointer::Val
        type(ResValinfo),pointer::ValHead,PVal,ValLast
        type(Resinfo),pointer::next
    end type
    
    type MeshGroupInfo
        character(70)   ::Name
        type(Groupinfo),dimension(:),pointer::Group
        type(Resinfo),dimension(:),pointer::Res
    endtype

    integer::mshunit,resunit,inpunit
    integer::ngroup,ncoor,nelem,nres,ntres,nstep,ndim,nzone
    logical::isMeshGroup,isOldFormat
    character(len=128)::arg
    integer::narg,iarg
    real::pCenter(3)
    real::angular
    character(len=128)::resType
    character(len=128)::fpath,text
    logical::eof
    
    
    
    type(Groupinfo),pointer::GroupHead,PGroup,GroupLast
    type(Groupinfo),dimension(:),pointer::Group
    type(coorinfo),pointer::CoorHead,PCoor,CoorLast
    type(coorinfo),dimension(:),pointer::Coor
    type(eleminfo),pointer::ElemHead,PElem,ElemLast
    type(eleminfo),dimension(:),pointer::Elem
    type(Resinfo),pointer::ResHead,PRes,ResLast
    type(Resinfo),dimension(:),pointer::Res
    type(Resinfo),dimension(:),pointer::ResTrans

    integer             ::ierr
    
    
    real,parameter:: pi=3.1415927

    end module
