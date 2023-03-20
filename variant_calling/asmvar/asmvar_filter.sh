#!/bin/bash

# Setting the current directory
cd variant_calling/asmvar/

# Increasing the number of files that can be opened at once
ulimit -n 2048

# Creating a variable for the reference genome
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# --- 1: Gathering and filtering all the vcf files for a given sample into files following pattern ${sample}_filtered.vcf

# DEPENDENCY: variant_calling/asmvar/ASMVAR_CALLING

# Creating a variable to hold the names of the 20 chromosomes
chrlist=$(seq -w 1 20 | xargs -I {} echo -n "Gm{} ")

# Creating a variable for the IDs of all the cultivars
# DEPENDENCY : ../../../srr_id_correspondence.txt
lines=$(cut -d " " -f1 ../../utilities/srr_id_correspondence.txt)

# The i variable iterates over all cultivars
for i in $lines
do
    # The j variable iterates over all chromosomes
    for j in $chrlist
    do
	    # DEPENDENCY : VCF files output by AsmVar
	    # Creating a variable for the name of the file
	    ij_file=${i}/asmvar_results_${j}.vcf
	    # Adding a header with all contig names based on the reference, and saving the result to a new file
	    # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
	    bcftools reheader --fai ${refgenome}.fai $ij_file > ${i}_${j}_header.vcf
    done
  
    # Concatenating the files for all chromosomes together, removing records without a "." FILTER tag, and removing SNPs
    bcftools concat -Ou $(ls ${i}_Gm??_header.vcf) | \
	    bcftools view --apply-filters .,PASS -Ou - | \
	    bcftools view -Ov --exclude-types snps - > ${i}_filtered.vcf

    # Deleting the temporary files with the added headers
    rm ${i}_Gm??_header.vcf
done

# --- 2: Normalizing the allelic representation in the VCF files of each sample

# Looping over the filtered vcf files for all lines
for i in $lines
do
	# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	bcftools norm -Ov -f $refgenome ${i}_filtered.vcf > ${i}_norm.vcf
done

# --- 3: Combining the variants from all samples into a single vcf file

bayestypertools combine -o asmvar_combined -v ESS:ESS_norm.vcf,HN001:HN001_norm.vcf,HN002:HN002_norm.vcf,HN003:HN003_norm.vcf,HN004:HN004_norm.vcf,HN005:HN005_norm.vcf,HN006:HN006_norm.vcf,HN007:HN007_norm.vcf,HN008:HN008_norm.vcf,HN009:HN009_norm.vcf,HN010:HN010_norm.vcf,HN012:HN012_norm.vcf,HN015:HN015_norm.vcf,HN016B:HN016B_norm.vcf,HN017B:HN017B_norm.vcf,HN018:HN018_norm.vcf,HN019:HN019_norm.vcf,HN021:HN021_norm.vcf,HN022:HN022_norm.vcf,HN023:HN023_norm.vcf,HN024:HN024_norm.vcf,HN025:HN025_norm.vcf,HN026:HN026_norm.vcf,HN028:HN028_norm.vcf,HN029:HN029_norm.vcf,HN030:HN030_norm.vcf,HN033:HN033_norm.vcf,HN035:HN035_norm.vcf,HN036:HN036_norm.vcf,HN038:HN038_norm.vcf,HN041:HN041_norm.vcf,HN042:HN042_norm.vcf,HN043:HN043_norm.vcf,HN044:HN044_norm.vcf,HN047:HN047_norm.vcf,HN048:HN048_norm.vcf,HN049:HN049_norm.vcf,HN051:HN051_norm.vcf,HN055:HN055_norm.vcf,HN056:HN056_norm.vcf,HN057:HN057_norm.vcf,HN058:HN058_norm.vcf,HN060:HN060_norm.vcf,HN061:HN061_norm.vcf,HN062:HN062_norm.vcf,HN064:HN064_norm.vcf,HN065:HN065_norm.vcf,HN066:HN066_norm.vcf,HN067:HN067_norm.vcf,HN068:HN068_norm.vcf,HN069:HN069_norm.vcf,HN070:HN070_norm.vcf,HN072:HN072_norm.vcf,HN073:HN073_norm.vcf,HN074:HN074_norm.vcf,HN075:HN075_norm.vcf,HN076:HN076_norm.vcf,HN077:HN077_norm.vcf,HN078:HN078_norm.vcf,HN079:HN079_norm.vcf,HN080:HN080_norm.vcf,HN082:HN082_norm.vcf,HN083:HN083_norm.vcf,HN084:HN084_norm.vcf,HN085:HN085_norm.vcf,HN086:HN086_norm.vcf,HN087:HN087_norm.vcf,HN088:HN088_norm.vcf,HN089:HN089_norm.vcf,HN091:HN091_norm.vcf,HN092:HN092_norm.vcf,HN093:HN093_norm.vcf,HN094:HN094_norm.vcf,HN095:HN095_norm.vcf,HN096:HN096_norm.vcf,HN097:HN097_norm.vcf,HN098:HN098_norm.vcf,HN099:HN099_norm.vcf,HN106:HN106_norm.vcf,HN107:HN107_norm.vcf,HNSM-25:HNSM-25_norm.vcf,HNY-29:HNY-29_norm.vcf,SRR1533279:SRR1533279_norm.vcf,SRR1533297:SRR1533297_norm.vcf,SRR1533345:SRR1533345_norm.vcf,SRR1533351:SRR1533351_norm.vcf,SRR1533362:SRR1533362_norm.vcf,SRR1533371:SRR1533371_norm.vcf,SRR1533384:SRR1533384_norm.vcf,SRR1533389:SRR1533389_norm.vcf,SRR1533395:SRR1533395_norm.vcf,SRR1533415:SRR1533415_norm.vcf,USB-002:USB-002_norm.vcf,USB-003:USB-003_norm.vcf,USB-013:USB-013_norm.vcf,USB-016:USB-016_norm.vcf,USB-017:USB-017_norm.vcf,USB-018:USB-018_norm.vcf,USB-020:USB-020_norm.vcf,USB-022:USB-022_norm.vcf,USB-023:USB-023_norm.vcf,USB-025:USB-025_norm.vcf,USB-026:USB-026_norm.vcf,USB-029:USB-029_norm.vcf,USB-032:USB-032_norm.vcf,USB-037:USB-037_norm.vcf,USB-038:USB-038_norm.vcf,USB-040:USB-040_norm.vcf,USB-041:USB-041_norm.vcf,USB-042:USB-042_norm.vcf,USB-044:USB-044_norm.vcf,USB-048:USB-048_norm.vcf,USB-050:USB-050_norm.vcf,USB-052:USB-052_norm.vcf,USB-059:USB-059_norm.vcf,USB-060:USB-060_norm.vcf,USB-061:USB-061_norm.vcf,USB-063:USB-063_norm.vcf,USB-064:USB-064_norm.vcf,USB-065:USB-065_norm.vcf,USB-066:USB-066_norm.vcf,USB-067:USB-067_norm.vcf,USB-068:USB-068_norm.vcf,USB-069:USB-069_norm.vcf,USB-070:USB-070_norm.vcf,USB-071:USB-071_norm.vcf,USB-072:USB-072_norm.vcf,USB-073:USB-073_norm.vcf,USB-074:USB-074_norm.vcf,USB-075:USB-075_norm.vcf,USB-076:USB-076_norm.vcf,USB-077:USB-077_norm.vcf,USB-079:USB-079_norm.vcf,USB-080:USB-080_norm.vcf,USB-081:USB-081_norm.vcf,USB-085:USB-085_norm.vcf,USB-086:USB-086_norm.vcf,USB-091:USB-091_norm.vcf,USB-092:USB-092_norm.vcf,USB-094:USB-094_norm.vcf,USB-095:USB-095_norm.vcf,USB-096:USB-096_norm.vcf,USB-097:USB-097_norm.vcf,USB-100:USB-100_norm.vcf,USB-102:USB-102_norm.vcf,USB-109:USB-109_norm.vcf,USB-116:USB-116_norm.vcf,USB-119:USB-119_norm.vcf,USB-127:USB-127_norm.vcf,USB-128:USB-128_norm.vcf,USB-143:USB-143_norm.vcf,USB-146:USB-146_norm.vcf,USB-147:USB-147_norm.vcf,USB-151:USB-151_norm.vcf,USB-153:USB-153_norm.vcf,USB-156:USB-156_norm.vcf,USB-158:USB-158_norm.vcf,USB-161:USB-161_norm.vcf,USB-162:USB-162_norm.vcf,USB-163:USB-163_norm.vcf,USB-167:USB-167_norm.vcf,USB-169:USB-169_norm.vcf,USB-171:USB-171_norm.vcf,USB-183:USB-183_norm.vcf,USB-185:USB-185_norm.vcf,USB-196:USB-196_norm.vcf,USB-206:USB-206_norm.vcf,USB-207:USB-207_norm.vcf,USB-208:USB-208_norm.vcf,USB-212:USB-212_norm.vcf,USB-214:USB-214_norm.vcf,USB-217:USB-217_norm.vcf,USB-221:USB-221_norm.vcf,USB-222:USB-222_norm.vcf,USB-223:USB-223_norm.vcf,USB-228:USB-228_norm.vcf,USB-233:USB-233_norm.vcf,USB-235:USB-235_norm.vcf,USB-236:USB-236_norm.vcf,USB-242:USB-242_norm.vcf,USB-245:USB-245_norm.vcf,USB-249:USB-249_norm.vcf,USB-250:USB-250_norm.vcf,USB-251:USB-251_norm.vcf,USB-252:USB-252_norm.vcf,USB-253:USB-253_norm.vcf,USB-257:USB-257_norm.vcf,USB-258:USB-258_norm.vcf,USB-259:USB-259_norm.vcf,USB-261:USB-261_norm.vcf,USB-271:USB-271_norm.vcf,USB-272:USB-272_norm.vcf,USB-274:USB-274_norm.vcf,USB-277:USB-277_norm.vcf,USB-282:USB-282_norm.vcf,USB-284:USB-284_norm.vcf,USB-285:USB-285_norm.vcf,USB-287:USB-287_norm.vcf,USB-298:USB-298_norm.vcf,USB-301:USB-301_norm.vcf,USB-302:USB-302_norm.vcf,USB-303:USB-303_norm.vcf,USB-304:USB-304_norm.vcf,USB-306:USB-306_norm.vcf,USB-307:USB-307_norm.vcf,USB-310:USB-310_norm.vcf,USB-312:USB-312_norm.vcf,USB-313:USB-313_norm.vcf,USB-314:USB-314_norm.vcf,USB-315:USB-315_norm.vcf,USB-316:USB-316_norm.vcf,USB-318:USB-318_norm.vcf,USB-319:USB-319_norm.vcf,USB-320:USB-320_norm.vcf,USB-322:USB-322_norm.vcf,USB-323:USB-323_norm.vcf,USB-324:USB-324_norm.vcf,USB-325:USB-325_norm.vcf,USB-326:USB-326_norm.vcf,USB-327:USB-327_norm.vcf,USB-330:USB-330_norm.vcf,USB-338:USB-338_norm.vcf,USB-341:USB-341_norm.vcf,USB-342:USB-342_norm.vcf,USB-345:USB-345_norm.vcf,USB-346:USB-346_norm.vcf,USB-347:USB-347_norm.vcf,USB-354:USB-354_norm.vcf,USB-357:USB-357_norm.vcf,USB-358:USB-358_norm.vcf,USB-363:USB-363_norm.vcf,USB-364:USB-364_norm.vcf,USB-366:USB-366_norm.vcf,USB-367:USB-367_norm.vcf,USB-368:USB-368_norm.vcf,USB-372:USB-372_norm.vcf,USB-375:USB-375_norm.vcf,USB-376:USB-376_norm.vcf,USB-377:USB-377_norm.vcf,USB-378:USB-378_norm.vcf,USB-381:USB-381_norm.vcf,USB-382:USB-382_norm.vcf,USB-383:USB-383_norm.vcf,USB-384:USB-384_norm.vcf,USB-385:USB-385_norm.vcf,USB-386:USB-386_norm.vcf,USB-390:USB-390_norm.vcf,USB-400:USB-400_norm.vcf,USB-412:USB-412_norm.vcf,USB-414:USB-414_norm.vcf,USB-418:USB-418_norm.vcf,USB-419:USB-419_norm.vcf,USB-421:USB-421_norm.vcf,USB-423:USB-423_norm.vcf,USB-424:USB-424_norm.vcf,USB-425:USB-425_norm.vcf,USB-430:USB-430_norm.vcf,USB-431:USB-431_norm.vcf,USB-432:USB-432_norm.vcf,USB-435:USB-435_norm.vcf,USB-436:USB-436_norm.vcf,USB-437:USB-437_norm.vcf,USB-438:USB-438_norm.vcf,USB-440:USB-440_norm.vcf,USB-441:USB-441_norm.vcf,USB-443:USB-443_norm.vcf,USB-444:USB-444_norm.vcf,USB-446:USB-446_norm.vcf,USB-447:USB-447_norm.vcf,USB-448:USB-448_norm.vcf,USB-449:USB-449_norm.vcf,USB-450:USB-450_norm.vcf,USB-451:USB-451_norm.vcf,USB-452:USB-452_norm.vcf,USB-453:USB-453_norm.vcf,USB-455:USB-455_norm.vcf,USB-459:USB-459_norm.vcf,USB-460:USB-460_norm.vcf,USB-462:USB-462_norm.vcf,USB-463:USB-463_norm.vcf,USB-465:USB-465_norm.vcf,USB-467:USB-467_norm.vcf,USB-468:USB-468_norm.vcf,USB-470:USB-470_norm.vcf,USB-472:USB-472_norm.vcf,USB-473:USB-473_norm.vcf,USB-475:USB-475_norm.vcf,USB-476:USB-476_norm.vcf,USB-477:USB-477_norm.vcf,USB-479:USB-479_norm.vcf,USB-480:USB-480_norm.vcf,USB-481:USB-481_norm.vcf,USB-484:USB-484_norm.vcf,USB-486:USB-486_norm.vcf,USB-497:USB-497_norm.vcf,USB-498:USB-498_norm.vcf,USB-500:USB-500_norm.vcf,USB-502:USB-502_norm.vcf,USB-504:USB-504_norm.vcf,USB-505:USB-505_norm.vcf,USB-512:USB-512_norm.vcf,USB-520:USB-520_norm.vcf,USB-524:USB-524_norm.vcf,USB-525:USB-525_norm.vcf,USB-526:USB-526_norm.vcf,USB-527:USB-527_norm.vcf,USB-528:USB-528_norm.vcf,USB-530:USB-530_norm.vcf,USB-531:USB-531_norm.vcf,USB-533:USB-533_norm.vcf,USB-534:USB-534_norm.vcf,USB-538:USB-538_norm.vcf,USB-541:USB-541_norm.vcf,USB-553:USB-553_norm.vcf,USB-560:USB-560_norm.vcf,USB-562:USB-562_norm.vcf,USB-565:USB-565_norm.vcf,USB-580:USB-580_norm.vcf,USB-584:USB-584_norm.vcf,USB-585:USB-585_norm.vcf,USB-593:USB-593_norm.vcf,USB-596:USB-596_norm.vcf,USB-597:USB-597_norm.vcf,USB-602:USB-602_norm.vcf,USB-604:USB-604_norm.vcf,USB-609:USB-609_norm.vcf,USB-613:USB-613_norm.vcf,USB-615:USB-615_norm.vcf,USB-616:USB-616_norm.vcf,USB-618:USB-618_norm.vcf,USB-619:USB-619_norm.vcf,USB-623:USB-623_norm.vcf,USB-630:USB-630_norm.vcf,USB-637:USB-637_norm.vcf,USB-642:USB-642_norm.vcf,USB-650:USB-650_norm.vcf,USB-656:USB-656_norm.vcf,USB-663:USB-663_norm.vcf,USB-664:USB-664_norm.vcf,USB-666:USB-666_norm.vcf,USB-669:USB-669_norm.vcf,USB-670:USB-670_norm.vcf,USB-672:USB-672_norm.vcf,USB-673:USB-673_norm.vcf,USB-674:USB-674_norm.vcf,USB-676:USB-676_norm.vcf,USB-677:USB-677_norm.vcf,USB-678:USB-678_norm.vcf,USB-679:USB-679_norm.vcf,USB-684:USB-684_norm.vcf,USB-688:USB-688_norm.vcf,USB-689:USB-689_norm.vcf,USB-694:USB-694_norm.vcf,USB-697:USB-697_norm.vcf,USB-699:USB-699_norm.vcf,USB-700:USB-700_norm.vcf,USB-701:USB-701_norm.vcf,USB-703:USB-703_norm.vcf,USB-711:USB-711_norm.vcf,USB-713:USB-713_norm.vcf,USB-729:USB-729_norm.vcf,USB-730:USB-730_norm.vcf,USB-737:USB-737_norm.vcf,USB-740:USB-740_norm.vcf,USB-745:USB-745_norm.vcf,USB-758:USB-758_norm.vcf,USB-759:USB-759_norm.vcf,USB-762:USB-762_norm.vcf,USB-763:USB-763_norm.vcf,USB-769:USB-769_norm.vcf,USB-773:USB-773_norm.vcf,USB-777:USB-777_norm.vcf,USB-778:USB-778_norm.vcf,USB-779:USB-779_norm.vcf,USB-780:USB-780_norm.vcf,USB-781:USB-781_norm.vcf,USB-782:USB-782_norm.vcf,USB-783:USB-783_norm.vcf,USB-784:USB-784_norm.vcf,USB-786:USB-786_norm.vcf,USB-792:USB-792_norm.vcf,USB-793:USB-793_norm.vcf,USB-794:USB-794_norm.vcf,USB-798:USB-798_norm.vcf,USB-801:USB-801_norm.vcf,USB-802:USB-802_norm.vcf,USB-803:USB-803_norm.vcf,USB-812:USB-812_norm.vcf,USB-819:USB-819_norm.vcf

# --- 4: Splitting the multiallelic records introduced by bayesTyperTools combine
bcftools norm --multiallelics -any -Ov asmvar_combined.vcf > asmvar_split.vcf

# --- 5: Adding the SVTYPE annotation along with the appropriate header line
# DEPENDENCY : scripts/add_svtype.awk
# DEPENDENCY : svtype_header_line.txt
../../scripts/add_svtype.awk asmvar_split.vcf | \
	awk 'BEGIN {OFS = "\t"} /^#/ {print} !/^#/ {sub("ACO=.*;", "ACO=asmvar;"); print}' | \
	bcftools annotate --header-lines svtype_header_line.txt -Ov - > asmvar_annotated.vcf

# --- 6: Extracting the SVs >= 50 bp from the annotated VCF file and setting the ID to the name of the caller + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk asmvar_annotated.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > asmvar_filtered.vcf

