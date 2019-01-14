clear
clc
load dataSet3.mat


warning('off');

        lncSim=miRNASS( LD_adjmat, disSim );

        disSim01  = GSD( LD_adjmat );
        lncSim01  = GSM( LD_adjmat );
        
        disSim02  = combineSim(disSim,disSim01);
        lncSim02  = combineSim(lncSim,lncSim01);
        
        KK=10;
        r=0.4;
        ld_adjmat_new=WKNKN( LD_adjmat, lncSim, disSim, KK, r );
        
        matPredict=NCPLDA(lncSim02, disSim02, ld_adjmat_new);


[NCP_rank,NCP_rank_known] =Rank_miRNAs( matPredict, LD_adjmat, lncRNA_Name, disease_Name);

Write_file( NCP_rank )
