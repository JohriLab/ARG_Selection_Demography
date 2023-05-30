grep -v MULTIALLELIC test.vcf >> bial_test.vcf
../../relate_v1.1.9_x86_64_dynamic/bin/RelateFileFormats --mode ConvertFromVcf --haps test.haps --sample test.sample -i test
#OR
../../relate_v1.1.9_x86_64_dynamic/bin/RelateFileFormats --mode ConvertFromVcf --haps bial_test.haps --sample bial_test.sample -i bial_test

../../relate_v1.1.9_MacOSX_M1/bin/Relate --mode All -m 1.25e-6 -N 2000 --haps bial_test.haps --sample bial_test.sample --map test.map --seed 1 -o bial_test

../../relate_v1.1.9_MacOSX_M1/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i bial_test \
              -m 1.25e-6 \
              --poplabels bial_test.poplabels \
              --seed 1 \
              -o bial_test_popsize
             
#bial_test_popsize.pdf has pre-made plot
#bial_test_popsize.coal can be converted with y=0.5/[line 3, removing first two values] ; x = [line 2]
