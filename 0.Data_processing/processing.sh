# list of gorilla individuals
fil="/scratch/devel/mkuhlwilm/arch/gorilla_noT.lst"
spec="gorilla"
/apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -h /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz > /scratch/devel/mkuhlwilm/arch/subsets/gehdr.txt
# for each chromosome $chrom, perform filtering
zcat /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz | sed -e 's/;MLEAC=.*;MLEAF=.*;/;/g' | /apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -a -U -V indels - | /apps/BCFTOOLS/1.9/bin/bcftools view -m 2 - | /apps/BCFTOOLS/1.9/bin/bcftools filter --exclude 'GT~"\."' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '(MIN(FMT/DP) >= 6) && (MAX(FMT/DP) <= 100) && (MIN(FMT/MQ) > 20)' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '((FMT/MQ0)/(FMT/DP)) < 0.1' | intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed | intersectBed -wa -v -header -a - -b <(zcat /scratch/devel/mkuhlwilm/arch/filter/Na_"$spec"_"$chrom".bed.gz ) | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | /apps/BCFTOOLS/1.9/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n' - | sed 's/\//|/g' | cat /scratch/devel/mkuhlwilm/arch/subsets/gehdr.txt - | bgzip -f > /scratch/devel/mkuhlwilm/arch/subsets/3_"$spec"_"$chrom".vcf.gz
