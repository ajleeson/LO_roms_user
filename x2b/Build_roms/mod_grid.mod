  �2  P   k820309    %          2021.1      -�^d                                                                                                          
       mod_grid.f90 MOD_GRID          @       �                                  
                                                                                                             #         @                                                       #NG    #LBI    #UBI    #LBJ    #UBJ    #LBIJ 	   #UBIJ 
                                                           
                                                       
                                                       
                                                       
                                                       
                                                       
                                  	                     
                                  
           #         @                                                       #NG              
                                             #         @                                                       #NG    #TILE    #MODEL                                                                
                                                       
                                                       
                                                               @                                '�             =      #ANGLER    #COSANGLER    #SINANGLER    #F    #FOMN    #GRDSCL    #H    #LATP    #LATR    #LATU    #LATV    #LONP    #LONR    #LONU    #LONV     #MYLON !   #OMN "   #OM_P #   #OM_R $   #OM_U %   #OM_V &   #ON_P '   #ON_R (   #ON_U )   #ON_V *   #PM +   #PN ,   #PMON_P -   #PMON_R .   #PMON_U /   #PMON_V 0   #PNOM_P 1   #PNOM_R 2   #PNOM_U 3   #PNOM_V 4   #ZOBOT 5   #RDRAG2 6   #XP 7   #XR 8   #XU 9   #XV :   #YP ;   #YR <   #YU =   #YV >   #HZ ?   #HUON @   #HVOM A   #Z0_R B   #Z0_W C   #Z_R D   #Z_V E   #Z_W F   #PMASK G   #RMASK H   #UMASK I   #VMASK J   #PMASK_FULL K   #RMASK_FULL L   #UMASK_FULL M   #VMASK_FULL N               �                                                           
            &                   &                                                       �                                         `                 
            &                   &                                                       �                                         �                 
            &                   &                                                       �                                                          
            &                   &                                                       �                                         �                
            &                   &                                                       �                                         �                
            &                   &                                                       �                                         @                
            &                   &                                                       �                                         �                
            &                   &                                                       �                                                       	   
            &                   &                                                       �                                         `             
   
            &                   &                                                       �                                         �                
            &                   &                                                       �                                                          
            &                   &                                                       �                                         �                
            &                   &                                                       �                                         �                
            &                   &                                                       �                                          @                
            &                   &                                                       �                             !            �                
            &                   &                                                       �                             "                             
            &                   &                                                       �                             #            `                
            &                   &                                                       �                             $            �                
            &                   &                                                       �                             %                             
            &                   &                                                       �                             &            �                
            &                   &                                                       �                             '            �                
            &                   &                                                       �                             (            @                
            &                   &                                                       �                             )            �                
            &                   &                                                       �                             *             	                
            &                   &                                                       �                             +            `	                
            &                   &                                                       �                             ,            �	                
            &                   &                                                       �                             -             
                
            &                   &                                                       �                             .            �
                
            &                   &                                                       �                             /            �
                
            &                   &                                                       �                             0            @                
            &                   &                                                       �                             1            �                 
            &                   &                                                       �                             2                          !   
            &                   &                                                       �                             3            `             "   
            &                   &                                                       �                             4            �             #   
            &                   &                                                       �                             5                          $   
            &                   &                                                       �                             6            �             %   
            &                   &                                                       �                             7            �             &   
            &                   &                                                       �                             8            @             '   
            &                   &                                                       �                             9            �             (   
            &                   &                                                       �                             :                          )   
            &                   &                                                       �                             ;            `             *   
            &                   &                                                       �                             <            �             +   
            &                   &                                                       �                             =                          ,   
            &                   &                                                       �                             >            �             -   
            &                   &                                                       �                             ?            �             .   
            &                   &                   &                                                       �                             @            X             /   
            &                   &                   &                                                       �                             A            �             0   
            &                   &                   &                                                       �                             B            H             1   
            &                   &                   &                                                       �                             C            �             2   
            &                   &                   &                                                       �                             D            8             3   
            &                   &                   &                                                       �                             E            �             4   
            &                   &                   &                                                       �                             F            (             5   
            &                   &                   &                                                       �                             G            �             6   
            &                   &                                                       �                             H                          7   
            &                   &                                                       �                             I            `             8   
            &                   &                                                       �                             J            �             9   
            &                   &                                                       �                             K                          :   
            &                   &                                                       �                             L            �             ;   
            &                   &                                                       �                             M            �             <   
            &                   &                                                       �                             N            @             =   
            &                   &                                                    @ @                               O            �                       &                                           #T_GRID       �         fn#fn    �   @   J   MOD_KINDS    �   p       R8+MOD_KINDS    n  �       ALLOCATE_GRID !   $  @   a   ALLOCATE_GRID%NG "   d  @   a   ALLOCATE_GRID%LBI "   �  @   a   ALLOCATE_GRID%UBI "   �  @   a   ALLOCATE_GRID%LBJ "   $  @   a   ALLOCATE_GRID%UBJ #   d  @   a   ALLOCATE_GRID%LBIJ #   �  @   a   ALLOCATE_GRID%UBIJ     �  P       DEALLOCATE_GRID #   4  @   a   DEALLOCATE_GRID%NG     t  �       INITIALIZE_GRID #     @   a   INITIALIZE_GRID%NG %   K  @   a   INITIALIZE_GRID%TILE &   �  @   a   INITIALIZE_GRID%MODEL    �  �      T_GRID    �  �   a   T_GRID%ANGLER !   G	  �   a   T_GRID%COSANGLER !   �	  �   a   T_GRID%SINANGLER    �
  �   a   T_GRID%F    K  �   a   T_GRID%FOMN    �  �   a   T_GRID%GRDSCL    �  �   a   T_GRID%H    O  �   a   T_GRID%LATP    �  �   a   T_GRID%LATR    �  �   a   T_GRID%LATU    S  �   a   T_GRID%LATV    �  �   a   T_GRID%LONP    �  �   a   T_GRID%LONR    W  �   a   T_GRID%LONU      �   a   T_GRID%LONV    �  �   a   T_GRID%MYLON    [  �   a   T_GRID%OMN      �   a   T_GRID%OM_P    �  �   a   T_GRID%OM_R    _  �   a   T_GRID%OM_U      �   a   T_GRID%OM_V    �  �   a   T_GRID%ON_P    c  �   a   T_GRID%ON_R      �   a   T_GRID%ON_U    �  �   a   T_GRID%ON_V    g  �   a   T_GRID%PM      �   a   T_GRID%PN    �  �   a   T_GRID%PMON_P    k  �   a   T_GRID%PMON_R      �   a   T_GRID%PMON_U    �  �   a   T_GRID%PMON_V    o  �   a   T_GRID%PNOM_P      �   a   T_GRID%PNOM_R    �  �   a   T_GRID%PNOM_U    s  �   a   T_GRID%PNOM_V       �   a   T_GRID%ZOBOT    �   �   a   T_GRID%RDRAG2    w!  �   a   T_GRID%XP    #"  �   a   T_GRID%XR    �"  �   a   T_GRID%XU    {#  �   a   T_GRID%XV    '$  �   a   T_GRID%YP    �$  �   a   T_GRID%YR    %  �   a   T_GRID%YU    +&  �   a   T_GRID%YV    �&  �   a   T_GRID%HZ    �'  �   a   T_GRID%HUON    _(  �   a   T_GRID%HVOM    #)  �   a   T_GRID%Z0_R    �)  �   a   T_GRID%Z0_W    �*  �   a   T_GRID%Z_R    o+  �   a   T_GRID%Z_V    3,  �   a   T_GRID%Z_W    �,  �   a   T_GRID%PMASK    �-  �   a   T_GRID%RMASK    O.  �   a   T_GRID%UMASK    �.  �   a   T_GRID%VMASK "   �/  �   a   T_GRID%PMASK_FULL "   S0  �   a   T_GRID%RMASK_FULL "   �0  �   a   T_GRID%UMASK_FULL "   �1  �   a   T_GRID%VMASK_FULL    W2  �       GRID 