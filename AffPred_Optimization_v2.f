      program main
      implicit none

      integer itime,itime_max
      parameter (itime_max=10000)
      integer ndata
      parameter (ndata=3205)
      integer ndata2
      parameter (ndata2=158)
      integer restype
      parameter (restype=20)
      integer CGAA_num
      parameter (CGAA_num=20)
      integer combination_num
      parameter (combination_num=1)
      real*8 NonInf_W,Inf_W
      parameter (NonInf_W=1.0,Inf_W=0.7)

      character*3 res_type(restype)
      integer res_ind(restype)
      integer Interface_Contact_Map(restype,restype,ndata)
      integer Interface_Contact_Map_HH(restype,restype,ndata)
      integer Interface_Contact_Map_SS(restype,restype,ndata)
      integer Interface_Contact_Map_LL(restype,restype,ndata)
      integer Interface_Contact_Map_HL(restype,restype,ndata)
      integer Interface_Contact_Map_SL(restype,restype,ndata)
      integer Interface_Contact_Map_HS(restype,restype,ndata)
      integer Noninterface_Res_H_S(restype,ndata)
      integer Noninterface_Res_S_S(restype,ndata)
      integer Noninterface_Res_L_S(restype,ndata)
      integer Noninterface_Res_H_C(restype,ndata)
      integer Noninterface_Res_S_C(restype,ndata)
      integer Noninterface_Res_L_C(restype,ndata)
      real*8 weight(restype,restype),weight_new(restype,restype)
      integer i,j,k,i2,j2,isample,i_comb
      real*8 calc_aff(ndata),norm_calc_aff(ndata)
      real*8 EXP_aff(ndata),train_exp_aff(ndata)
      real*8 norm_train_exp_aff(ndata)
      real*8 optimize_aff(ndata)
      real*8 test_exp_aff(ndata)
      character*4 PDBid(ndata)
      real*8 cc,cc_max
      integer pick_i,pick_j,pick_k,pick_mode
      integer train_Inf_sample(restype,restype,ndata)
      integer train_Inf_sample_HH(restype,restype,ndata)
      integer train_Inf_sample_SS(restype,restype,ndata)
      integer train_Inf_sample_LL(restype,restype,ndata)
      integer train_Inf_sample_HL(restype,restype,ndata)
      integer train_Inf_sample_SL(restype,restype,ndata)
      integer train_Inf_sample_HS(restype,restype,ndata)

      integer train_Noninf_sample_H_S(restype,ndata)
      integer train_Noninf_sample_S_S(restype,ndata)
      integer train_Noninf_sample_L_S(restype,ndata)
      integer train_Noninf_sample_H_C(restype,ndata)
      integer train_Noninf_sample_S_C(restype,ndata)
      integer train_Noninf_sample_L_C(restype,ndata)

      integer test_Inf_sample(ndata,restype,restype)
      integer test_Inf_sample_HH(ndata,restype,restype)
      integer test_Inf_sample_SS(ndata,restype,restype)
      integer test_Inf_sample_LL(ndata,restype,restype)
      integer test_Inf_sample_HL(ndata,restype,restype)
      integer test_Inf_sample_SL(ndata,restype,restype)
      integer test_Inf_sample_HS(ndata,restype,restype)

      integer test_Noninf_sample_H_S(ndata,restype)
      integer test_Noninf_sample_S_S(ndata,restype)
      integer test_Noninf_sample_L_S(ndata,restype)
      integer test_Noninf_sample_H_C(ndata,restype)
      integer test_Noninf_sample_S_C(ndata,restype)
      integer test_Noninf_sample_L_C(ndata,restype)

      real*8 test_aff(ndata)
      integer index
      integer CGAA_index(restype)
      real*8 weight_CG(CGAA_num,CGAA_num)
      real*8 weight_CG_new(CGAA_num,CGAA_num)
      real*8 weight_Inf(CGAA_num,CGAA_num,6)
      real*8 weight_Inf_new(CGAA_num,CGAA_num,6)
      real*8 weight_Noninf(CGAA_num,6)
      real*8 weight_Noninf_new(CGAA_num,6)

      real*8 ave_train_cc,train_cc(ndata2)
      integer res_type_flag(restype),type_index
      integer index_testnum,index_trainnum
      integer test_isample,train_isample
      character*1 IP
      character*4 test_ip(ndata),train_ip(ndata)
      character*3 train_id(ndata2)
      real*8 Ave,SD
      character*4 recna
      integer in
      integer train_flag

      real rand3
      double precision r3      
      real rand4
      double precision r4

      r3=5.0
      r4=5.0

ccccccccccccccccccccccccc

      open (unit=26,file='SKEMPI_training/'//
     &     'trainingset_list.dat',status='old')
      
      do i=1,ndata2
         read(26,6000) train_id(i)
      enddo
 6000 format(A3)
      close(26)
      
      open (unit=26,file='res_index.dat',status='old')
      do i=1,restype
         read(26,7000) res_type(i),res_ind(i)
      enddo
 7000 format(A3,I3)
      
      close(26)


      do i=1,ndata
         do i2=1,restype
            do j2=1,restype
               Interface_Contact_Map(i2,j2,i)=0
               Interface_Contact_Map_HH(i2,j2,i)=0
               Interface_Contact_Map_SS(i2,j2,i)=0
               Interface_Contact_Map_LL(i2,j2,i)=0
               Interface_Contact_Map_HL(i2,j2,i)=0
               Interface_Contact_Map_SL(i2,j2,i)=0
               Interface_Contact_Map_HS(i2,j2,i)=0
            enddo
         enddo
      enddo
      do i=1,ndata
         do i2=1,restype
            Noninterface_Res_H_S(i2,i)=0
            Noninterface_Res_S_S(i2,i)=0
            Noninterface_Res_L_S(i2,i)=0
            Noninterface_Res_H_C(i2,i)=0
            Noninterface_Res_S_C(i2,i)=0
            Noninterface_Res_L_C(i2,i)=0
         enddo
      enddo

      open (unit=10,file=
     &     'SKEMPI_InterfaceStat_03262018.dat',
     &     status='old')

      do i=1,ndata
         
         read(10,2200) PDBid(i),EXP_aff(i)
 2200    format(A4,1x,F6.2)
         
         do i2=1,restype
            do j2=1,i2
               read(10,2300) Interface_Contact_Map(i2,j2,i),
     &              Interface_Contact_Map_HH(i2,j2,i),
     &              Interface_Contact_Map_SS(i2,j2,i),
     &              Interface_Contact_Map_LL(i2,j2,i),
     &              Interface_Contact_Map_HL(i2,j2,i),
     &              Interface_Contact_Map_SL(i2,j2,i),
     &              Interface_Contact_Map_HS(i2,j2,i)
            enddo
         enddo
         do i2=1,restype
            read(10,2400) Noninterface_Res_H_S(i2,i),
     &           Noninterface_Res_S_S(i2,i),
     &           Noninterface_Res_L_S(i2,i),
     &           Noninterface_Res_H_C(i2,i),
     &           Noninterface_Res_S_C(i2,i),
     &           Noninterface_Res_L_C(i2,i)
         enddo
 2300    format(7x,7(1x,I10))
 2400    format(3x,6(1x,I10))
      enddo
      close(10)


c      do i=1,ndata
c         test_exp_aff(i)=0
c         do i2=1,restype
c            do j2=1,i2               
c               test_exp_aff(i)=test_exp_aff(i)+
c     &              Interface_Contact_Map(i2,j2,i)*
c     &              (-1.0)
c            enddo
c         enddo
c      enddo
c      cc=0
c      call CalCorrCoff(EXP_aff,test_exp_aff,ndata,cc)
c      print*,cc
    

      do i_comb=1,combination_num

         do i=1,restype
c            res_type_flag(i)=0
            CGAA_index(i)=i
         enddo

c         do i=1,CGAA_num
c 2002       continue
c            type_index=int(rand3(r3)*restype)+1
c            if(res_type_flag(type_index).eq.0)then
c               CGAA_index(type_index)=i
c               res_type_flag(type_index)=1
c            else
c               goto 2002
c            endif
c         enddo

c         do i=1,restype
c            if(res_type_flag(i).eq.0)then
c               CGAA_index(i)=int(rand3(r3)*CGAA_num)+1
c            endif
c         enddo

c         do i=1,restype
c            print*,i,CGAA_index(i)
c         enddo
         

c         CGAA_index(1)=int(rand3(r3)*CGAA_num)+1 !GLY
c         CGAA_index(2)=int(rand3(r3)*CGAA_num)+1 !ALA
c         CGAA_index(3)=int(rand3(r3)*CGAA_num)+1 !VAL
c         CGAA_index(4)=int(rand3(r3)*CGAA_num)+1 !ILE
c         CGAA_index(5)=int(rand3(r3)*CGAA_num)+1 !LEU
c         CGAA_index(6)=int(rand3(r3)*CGAA_num)+1 !SER
c         CGAA_index(7)=int(rand3(r3)*CGAA_num)+1 !THR
c         CGAA_index(8)=int(rand3(r3)*CGAA_num)+1 !ASP
c         CGAA_index(9)=int(rand3(r3)*CGAA_num)+1 !ASN
c         CGAA_index(10)=int(rand3(r3)*CGAA_num)+1 !GLU
c         CGAA_index(11)=int(rand3(r3)*CGAA_num)+1 !GLN
c         CGAA_index(12)=int(rand3(r3)*CGAA_num)+1 !LYS
c         CGAA_index(13)=int(rand3(r3)*CGAA_num)+1 !ARG
c         CGAA_index(14)=int(rand3(r3)*CGAA_num)+1 !CYS
c         CGAA_index(15)=int(rand3(r3)*CGAA_num)+1 !MET
c         CGAA_index(16)=int(rand3(r3)*CGAA_num)+1 !PHE
c         CGAA_index(17)=int(rand3(r3)*CGAA_num)+1 !TYR
c         CGAA_index(18)=int(rand3(r3)*CGAA_num)+1 !TRP
c         CGAA_index(19)=int(rand3(r3)*CGAA_num)+1 !HIS
c         CGAA_index(20)=int(rand3(r3)*CGAA_num)+1 !PRO

         index_testnum=0
         
         do isample=1,ndata2

            test_isample=0
            train_isample=0

            open (unit=1,file='SKEMPI_training/'//
     &           'train_'//train_id(isample)//'.list',status='old',
     &           access='append')
            
            write(1,800) 'ENDD'
 800        format(A4)
            
            close(1)

            open (unit=1,file='SKEMPI_training/'//
     &           'train_'//train_id(isample)//'.list',status='old',
     &           access='sequential')
            in=0

 900        read(1,5000) recna
 5000       format(A4)

            in=in+1

            if(recna.eq.'ENDD')then
               goto 910
            endif

            goto 900

 910        continue
            
            close(1)
            
            test_isample=in-2

            open (unit=26,file='SKEMPI_training/'//
     &           'train_'//train_id(isample)//'.list',status='old')
            read(26,*) 
            do i=1,test_isample
               read(26,8000) test_ip(i)
 8000          format(A4)
            enddo
            close(26)


            do i=1,test_isample
               do j=1,ndata
                  if(test_ip(i).eq.PDBid(j))then
                     index_testnum=index_testnum+1
                     test_exp_aff(index_testnum)=EXP_aff(j)
                     do i2=1,restype
                        do j2=1,i2
                           test_Inf_sample(index_testnum,i2,j2)=
     &                          Interface_Contact_Map(i2,j2,j)
                           test_Inf_sample_HH(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_HH(i2,j2,j)
                           test_Inf_sample_SS(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_SS(i2,j2,j)
                           test_Inf_sample_LL(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_LL(i2,j2,j)
                           test_Inf_sample_HL(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_HL(i2,j2,j)
                           test_Inf_sample_SL(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_SL(i2,j2,j)
                           test_Inf_sample_HS(index_testnum,i2,j2)=
     &                          Interface_Contact_Map_HS(i2,j2,j)
                        enddo
                     enddo
                     do i2=1,restype
                        test_Noninf_sample_H_S(index_testnum,i2)=
     &                       Noninterface_Res_H_S(i2,j)
                        test_Noninf_sample_S_S(index_testnum,i2)=
     &                       Noninterface_Res_S_S(i2,j)
                        test_Noninf_sample_L_S(index_testnum,i2)=
     &                       Noninterface_Res_L_S(i2,j)
                        test_Noninf_sample_H_C(index_testnum,i2)=
     &                       Noninterface_Res_H_C(i2,j)
                        test_Noninf_sample_S_C(index_testnum,i2)=
     &                       Noninterface_Res_S_C(i2,j)
                        test_Noninf_sample_L_C(index_testnum,i2)=
     &                       Noninterface_Res_L_C(i2,j)
                     enddo
                  endif
               enddo
            enddo

            index_trainnum=0
            do i=1,ndata
            	 train_flag=0		
               do j=1,test_isample
                  if(test_ip(j).eq.PDBid(i))then
                     train_flag=1
                  endif
               enddo
               if(train_flag.eq.0)then	
                  index_trainnum=index_trainnum+1
                  do i2=1,restype
                     do j2=1,i2
                        train_Inf_sample(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map(i2,j2,i)
                        train_Inf_sample_HH(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_HH(i2,j2,i)
                        train_Inf_sample_SS(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_SS(i2,j2,i)
                        train_Inf_sample_LL(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_LL(i2,j2,i)
                        train_Inf_sample_HL(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_HL(i2,j2,i)
                        train_Inf_sample_SL(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_SL(i2,j2,i)
                        train_Inf_sample_HS(i2,j2,index_trainnum)=
     &                       Interface_Contact_Map_HS(i2,j2,i)
                     enddo
                  enddo
                  do i2=1,restype
                     train_Noninf_sample_H_S(i2,index_trainnum)=
     &                    Noninterface_Res_H_S(i2,i)
                     train_Noninf_sample_S_S(i2,index_trainnum)=
     &                    Noninterface_Res_S_S(i2,i)
                     train_Noninf_sample_L_S(i2,index_trainnum)=
     &                    Noninterface_Res_L_S(i2,i)
                     train_Noninf_sample_H_C(i2,index_trainnum)=
     &                    Noninterface_Res_H_C(i2,i)
                     train_Noninf_sample_S_C(i2,index_trainnum)=
     &                    Noninterface_Res_S_C(i2,i)
                     train_Noninf_sample_L_C(i2,index_trainnum)=
     &                    Noninterface_Res_L_C(i2,i)
                  enddo
                  train_exp_aff(index_trainnum)=EXP_aff(i)
               endif
            enddo
            train_isample=index_trainnum
      

c            do i=1,train_isample
c               norm_train_exp_aff(i)=(train_exp_aff(i)-
c     &              Ave(train_exp_aff,train_isample))/
c     &              SD(train_exp_aff,train_isample)
c            enddo
            
            do i=1,ndata
               optimize_aff(i)=0
            enddo

            do i=1,restype
               do j=1,restype
                  do k=1,6
                     weight_Inf(i,j,k)=0
                  enddo
               enddo
            enddo
            do i=1,restype
               do k=1,6
                  weight_Noninf(i,j)=0
               enddo
            enddo
            
c            do i=1,CGAA_num
c               do j=1,CGAA_num
c                  weight_CG(i,j)=0
c               enddo
c            enddo
            
            cc_max=0.0
            
            do itime=1,itime_max
               
               do i=1,restype
                  do j=1,i
                     do k=1,6
                        weight_Inf_new(i,j,k)=weight_Inf(i,j,k)
                     enddo
                  enddo
               enddo
               do i=1,restype
                  do k=1,6
                     weight_Noninf_new(i,k)=weight_Noninf(i,k)
                  enddo
               enddo

               pick_mode=int(rand4(r4)*2.0)+1


               if(pick_mode.eq.1)then
                  pick_i=int(rand3(r3)*restype)+1
                  pick_j=int(rand4(r4)*pick_i)+1
                  pick_k=int(rand3(r3)*6.0)+1
                  
                  weight_Inf_new(pick_i,pick_j,pick_k)=Inf_W*
     &                 (rand4(r4)*2-1.0)
                  
               elseif(pick_mode.eq.2)then
                  
                  pick_i=int(rand3(r3)*restype)+1
                  pick_k=int(rand4(r4)*6.0)+1
                  
                  weight_Noninf_new(pick_i,pick_k)=NonInf_W*
     &                 (rand4(r4)*2-1.0)

               endif
               
               do i=1,train_isample
                  calc_aff(i)=0
                  do i2=1,restype
                     do j2=1,i2
                        calc_aff(i)=calc_aff(i)+
     &                       train_Inf_sample_HH(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,1)+
     &                       train_Inf_sample_SS(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,2)+
     &                       train_Inf_sample_LL(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,3)+
     &                       train_Inf_sample_HL(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,4)+
     &                       train_Inf_sample_SL(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,5)+
     &                       train_Inf_sample_HS(i2,j2,i)*
     &                       weight_Inf_new(i2,j2,6)                   
                     enddo
                  enddo
                  do i2=1,restype
                     calc_aff(i)=calc_aff(i)+
     &                    train_Noninf_sample_H_S(i2,i)*
     &                    weight_Noninf_new(i2,1)+
     &                    train_Noninf_sample_S_S(i2,i)*
     &                    weight_Noninf_new(i2,2)+
     &                    train_Noninf_sample_L_S(i2,i)*
     &                    weight_Noninf_new(i2,3)+
     &                    train_Noninf_sample_H_C(i2,i)*
     &                    weight_Noninf_new(i2,4)+
     &                    train_Noninf_sample_S_C(i2,i)*
     &                    weight_Noninf_new(i2,5)+
     &                    train_Noninf_sample_L_C(i2,i)*
     &                    weight_Noninf_new(i2,6)
                  enddo
               enddo
               
c               do i=1,train_isample
c                  norm_calc_aff(i)=(calc_aff(i)
c     &                 -Ave(calc_aff,train_isample))/
c     &                 SD(calc_aff,train_isample
c     &                 )
c               enddo

               cc=0
               call CalCorrCoff(train_exp_aff,calc_aff,
     &              train_isample,cc)
               
               
               
               if(cc.gt.cc_max)then
                  
                  cc_max=cc
                  
                  do i=1,restype
                     do j=1,i
                        do k=1,6
                           weight_Inf(i,j,k)=weight_Inf_new(i,j,k)
                        enddo
                     enddo
                     do k=1,6
                        weight_Noninf(i,k)=weight_Noninf_new(i,k)
                     enddo
                  enddo

                  do i=1,train_isample
                     optimize_aff(i)=calc_aff(i)
                  enddo
c                  print*,isample,itime,cc_max
                  
               endif
            
            enddo

c      do i=1,restype
c         do j=1,i
c            print*,res_type(i),res_type(j),weight(i,j)
c         enddo
c      enddo

c      do i=1,ndata
c            print*,i,EXP_aff(i),calc_aff(i)
c      enddo

c            do i=1,CGAA_num
c               do j=1,i
c                  print*,'weight',weight_CG(i,j)
c               enddo
c            enddo

            open (unit=10,file=
     &           'SKEMPI_AffPred_IW10NIW07_mc10000_WeightRec.dat'
     &           ,status='unknown'
     &           ,access='append')
            do i2=1,restype
               do j2=1,i2
                  write(10,1980) res_type(i2),res_type(j2),
     &                 weight_Inf(i2,j2,1),weight_Inf(i2,j2,2),
     &                 weight_Inf(i2,j2,3),weight_Inf(i2,j2,4),
     &                 weight_Inf(i2,j2,5),weight_Inf(i2,j2,6)
               enddo
            enddo
 1980       format(A3,1x,A3,1x,6(1x,F10.5))
            do i2=1,restype
               write(10,1981) res_type(i2),
     &              weight_Noninf(i2,1),weight_Noninf(i2,2),
     &              weight_Noninf(i2,3),weight_Noninf(i2,4),
     &              weight_Noninf(i2,5),weight_Noninf(i2,6)
            enddo
 1981       format(A3,1x,6(1x,F10.5))
            close(10)

            do i=1,test_isample
               test_aff(index_testnum-test_isample+i)=0
               do i2=1,restype
                  do j2=1,i2
                     test_aff(index_testnum-test_isample+i)=
     &                    test_aff(index_testnum-test_isample+i)+
     &                    test_Inf_sample_HH(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,1)+
     &                    test_Inf_sample_SS(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,2)+
     &                    test_Inf_sample_LL(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,3)+
     &                    test_Inf_sample_HL(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,4)+
     &                    test_Inf_sample_SL(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,5)+
     &                    test_Inf_sample_HS(index_testnum-
     &                    test_isample+i,i2,j2)*
     &                    weight_Inf(i2,j2,6)
                  enddo
               enddo
               do i2=1,restype
                  test_aff(index_testnum-test_isample+i)=
     &                 test_aff(index_testnum-test_isample+i)+
     &                 test_Noninf_sample_H_S(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,1)+
     &                 test_Noninf_sample_S_S(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,2)+
     &                 test_Noninf_sample_L_S(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,3)+
     &                 test_Noninf_sample_H_C(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,4)+
     &                 test_Noninf_sample_S_C(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,5)+
     &                 test_Noninf_sample_L_C(index_testnum-
     &                 test_isample+i,i2)*
     &                 weight_Noninf(i2,6)
               enddo
               test_aff(index_testnum-test_isample+i)=
     &              ((test_aff(index_testnum-test_isample+i)-
     &              Ave(optimize_aff,train_isample))/
     &              SD(optimize_aff,train_isample))*
     &              SD(train_exp_aff,train_isample)
     &              +Ave(train_exp_aff,train_isample)
               open (unit=10,file=
     &              'SKEMPI_AffPred_IW10NIW07_mc10000_v2.dat'
     &              ,status='unknown'
     &              ,access='append')
               write(10,1977) test_ip(i),
     &              test_exp_aff(index_testnum-test_isample+i),
     &              test_aff(index_testnum-test_isample+i)
 1977          format(A4,1x,F10.5,1x,F10.5)
               close(10)
            enddo
            
c           print*,'max cc ',cc_max

            train_cc(isample)=cc_max
            
c           print*,isample,test_exp_aff(isample),test_aff(isample)

         enddo

         ave_train_cc=0
         do isample=1,ndata2
            ave_train_cc=ave_train_cc+train_cc(isample)
         enddo
         ave_train_cc=ave_train_cc/real(ndata2)

         cc=0
         call CalCorrCoff(test_exp_aff,test_aff,ndata,cc)

         open (unit=10,file='SKEMPI_AffPred_IW10NIW07_mc10000_v2.dat'
     &        ,status='unknown'
     &        ,access='append')
c         do i=1,restype
c            write(10,*) i,' ',res_type(i),' ',CGAA_index(i)
c         enddo
         write(10,*) i_comb,' train cc ',ave_train_cc,' validate cc ',cc
         close(10)

      enddo

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccc

      subroutine CalCorrCoff(x,y,ndp,cc)
      implicit none

      integer ndp
      real*8 x(ndp)
      real*8 y(ndp)
      real*8 cc

      integer i
      real*8 ave_x,ave_y
      real*8 sd_x,sd_y
      real*8 PCC


c>>  calculate average of x
      
      ave_x=0
      do i=1,ndp
         ave_x=ave_x+x(i)
      enddo
      ave_x=ave_x/real(ndp)

c>>  calculate average of y
      
      ave_y=0
      do i=1,ndp
         ave_y=ave_y+y(i)
      enddo
      ave_y=ave_y/real(ndp)

c>>  calculate SD of x
      
      sd_x=0
      do i=1,ndp
         sd_x=sd_x+(x(i)-ave_x)**2
      enddo
      sd_x=sqrt(sd_x/real(ndp))

c>>  calculate SD of y
      
      sd_y=0
      do i=1,ndp
         sd_y=sd_y+(y(i)-ave_y)**2
      enddo
      sd_y=sqrt(sd_y/real(ndp))

c>>  calculate Pearson Correlation Coefficient

      PCC=0
      do i=1,ndp
         PCC=PCC+(x(i)-ave_x)*(y(i)-ave_y)
      enddo
      PCC=PCC/(sd_x*sd_y*real(ndp))

      cc=PCC

      return 
      end

c>>>>>>>>>>>>>>>>>>>>>>>>>.
c>>>>>>>>>>>>>>>>>>>>>>>>>.

      Function Ave(x,ndp)
      implicit none

      real*8 Ave
      integer ndp
      real*8 x(ndp)

      integer i
      real*8 ave_x

c>>  calculate average of x

      ave_x=0
      do i=1,ndp
         ave_x=ave_x+x(i)
      enddo
      ave_x=ave_x/real(ndp)

      Ave=ave_x

      end

cccccccccccccccccccccccccccccccccc

c>>>>>>>>>>>>>>>>>>>>>>>>>.

      Function SD(x,ndp)
      implicit none

      real*8 SD
      integer ndp
      real*8 x(ndp)

      integer i
      real*8 ave_x
      real*8 sd_x

c>>  calculate average of x
      
      ave_x=0
      do i=1,ndp
         ave_x=ave_x+x(i)
      enddo
      ave_x=ave_x/real(ndp)

c>>  calculate SD of x
      
      sd_x=0
      do i=1,ndp
         sd_x=sd_x+(x(i)-ave_x)**2
      enddo
      sd_x=sqrt(sd_x/real(ndp))

      SD=sd_x

      end

cccccccccccccccccccccccccccccccccc

