verbosity = 3
#: Verbosity level (0=errors, 1=warnings, 2=infos, 3=hints, debugs=4)

logfile = ""
#: Name of logfile. By default is set to '' and writes to standard output.

tools = ""
#: Path to the directory where external tools are. By default is set to '$PATH'.

refs = ["hg19", "hg38", "mm10", "m10x"]
HUMAN = "/fdb/igenomes/Homo_sapiens/UCSC"
MOUSE = "/fdb/igenomes/Mus_musculus/UCSC"
GATK = "/fdb/GATK_resource_bundle"

hg19 = {}
hg19["ref"] = f"{HUMAN}/hg19/Sequence/WholeGenomeFasta/genome.fa"
hg19["bwa_ref"] = f"{HUMAN}/hg19/Sequence/BWAIndex/genome.fa"
hg19["annot"] = f"{HUMAN}/hg19/Annotation/Genes/genes.gtf"
hg19["known_sites"] = [
    f"{GATK}/hg19/dbsnp_138.hg19.vcf.gz",
    f"{GATK}/hg19/1000G_phase1.snps.high_confidence.hg19.vcf.gz",
    f"{GATK}/hg19/1000G_phase1.indels.hg19.vcf.gz",
    f"{GATK}/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
]
hg19["defuse_config"] = "/usr/local/apps/defuse/0.8.1/config_hg19_ens69.txt"
hg19["defuse_db"] = "/fdb/defuse/hg19_ens69_newgmap"

hg38 = {}
hg38["ref"] = f"{HUMAN}/hg38/Sequence/WholeGenomeFasta/genome.fa"
hg38["bwa_ref"] = f"{HUMAN}/hg38/Sequence/BWAIndex/genome.fa"
hg38["annot"] = f"{HUMAN}/hg38/Annotation/Genes/genes.gtf"
hg38["known_sites"] = [
    f"{GATK}/hg38/dbsnp_138.hg38.vcf.gz",
    f"{GATK}/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    f"{GATK}/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
]

mm10 = {}
mm10["ref"] = f"{MOUSE}/mm10/Sequence/WholeGenomeFasta/genome.fa"
mm10["bwa_ref"] = f"{MOUSE}/mm10/Sequence/BWAIndex/genome.fa"
mm10["annot"] = f"{MOUSE}/mm10/Annotation/Genes/genes.gtf"
mm10["known_sites"] = [
    f"{GATK}/mm10/dbsnp146_fixedNames.vcf.gz",
    f"{GATK}/mm10/mgp.v5.merged.indels.dbSNP142.normed.fixedNames.vcf.gz",
    f"{GATK}/mm10/mgp.v5.merged.snps_all.dbSNP142.fixedNames.vcf.gz",
]
mm10["defuse_config"] = "/usr/local/apps/defuse/0.8.1/config_mm10_ens84.txt"
mm10["defuse_db"] = "/fdb/defuse/mm10_ens84_newgmap"

m10x = {}
m10x["ref"] = "/data/rashidimehrabf2/_database/refdata-gex-mm10/fasta/genome.fa"

tools = {}
tools["email"] = "farid.rsh@gmail.com"
tools["bwa"] = "bwa/0.7.17"
tools["samtools"] = "samtools/1.11"
tools["bamtools"] = "bamtools/2.5.1"
tools["bcftools"] = "bcftools/1.9"
tools["fastp"] = "fastp/0.20.1"
tools["star"] = "STAR/2.7.3a"
tools["gatk"] = "GATK/4.2.1.0"
tools["rsem"] = "rsem/1.3.2"
tools["snpeff"] = "snpEff/5.0"
tools["defuse"] = "defuse/0.8.1"
tools["bamreadcount"] = "bamreadcount/0.8.0"
tools["vartrix"] = "vartrix/1.1.22"
tools["strelka"] = "strelka/2.9.10"
tools["varscan"] = "varscan/2.4.3"
tools["sequenza"] = "sequenza-utils/3.0.0"
tools["pyclone"] = "pyclone/0.13.1"
