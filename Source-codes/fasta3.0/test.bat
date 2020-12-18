rem "starting fasta34 - protein"
fasta34 -q -S mgstm1.aa q > test_m1.ok2S
fasta34 -q -O test_m1p.ok2 -s P250 mgstm1.aa s
rem "done"
rem "starting fastx3"
fastx34 -q -S mgstm1.esq q > test_m1.ok2x
rem "done"
rem "starting fasty3"
fasty34 -q -S mgstm1.esq a > test_m1.ok2y
rem "done"
rem "starting fasta34 - DNA "
fasta34 -q mgstm1.seq %%MB 4 > test_m1.ok4
rem "done"
rem "starting ssearch3"
ssearch34 -q mgstm1.aa q > test_m1.ss 
ssearch34 -q -s P250 -O test_m1_p.ss mgstm1.aa s
rem "done"
rem "starting tfasta3"
tfasta34 -q mgstm1.aa %%MB > test_m1.tk2 
rem "done"
rem "starting tfastxy3"
tfastx34 -q -m 6 mgstm1.aa %%MB > test_m1_tx2.html
tfasty34 -q mgstm1.aa %%MB > test_m1.ty2 
rem "done"
