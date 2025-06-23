#!/bin/bash

ARRIBA_PATH="/home/bioinfo/tools/arriba_v2.4.0/run_arriba.sh"
REF_STAR="/home/bioinfo/RefSeq/RefSTAR/"
GTF_FILE="/home/bioinfo/RefSeq/Homo_sapiens.GRCh38.112.gtf"
FASTA_FILE="/home/bioinfo/RefSeq/RefBWA/hg38_v0_Homo_sapiens_assembly38.fasta"
BLACKLIST="/home/bioinfo/tools/arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv"
KNOWN_FUSIONS="/home/bioinfo/tools/arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv"
PROTEIN_DOMAINS="/home/bioinfo/tools/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"

DRAW_SCRIPT="/home/bioinfo/tools/arriba_v2.4.0/draw_fusions.R"
CYTOBANDS="/home/bioinfo/tools/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv"

FASTQ_DIR="$PWD"

for R1_FILE in ${FASTQ_DIR}/*R1_001.fastq.gz; do
  # Gera o nome do arquivo R2 correspondente
  R2_FILE="${R1_FILE/_R1_/_R2_}"

  # Verifica se o arquivo R2 correspondente existe
  if [[ -f "${R2_FILE}" ]]; then
    echo "🧬 Processando: ${R1_FILE} e ${R2_FILE}"

    # Executa o script Arriba
    ${ARRIBA_PATH} \
      ${REF_STAR} \
      ${GTF_FILE} \
      ${FASTA_FILE} \
      ${BLACKLIST} \
      ${KNOWN_FUSIONS} \
      ${PROTEIN_DOMAINS} 21 ${R1_FILE} ${R2_FILE}

    # Renomeia o arquivo fusions.tsv
    sample_name=$(basename "$R1_FILE" | cut -d'_' -f1)
    mv fusions.tsv "fusions_${sample_name}.tsv"

    # Rodar a visualização das fusões
    echo "🧬 Gerando figura para ${sample_name}..."
    ${DRAW_SCRIPT} \
      --fusions="fusions_${sample_name}.tsv" \
      --alignments=Aligned.sortedByCoord.out.bam \
      --output="fusions_${sample_name}.pdf" \
      --annotation=${GTF_FILE} \
      --cytobands=${CYTOBANDS} \
      --proteinDomains=${PROTEIN_DOMAINS}

    echo "✅ Finalizado ${sample_name}!"
  else
    echo "❌ Arquivo correspondente R2 não encontrado para: ${R1_FILE}"
  fi
done

echo "✅✅ Pipeline completo concluído para todos os pares!"
