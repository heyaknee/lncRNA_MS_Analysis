#!/bin/bash
SRA_DIR="path"  # Path to SRR files
OUTPUT_DIR="${SRA_DIR}/FASTQ_Output"  # Output directory for FASTQ files
MD5_FILE="${SRA_DIR}/md5_checksums_pairs.txt"  # MD5 checksum output file

# Create output directory if not exists
mkdir -p $OUTPUT_DIR

# List of remaining SRR files (single SRRs first)
SRR_FILES=(
    "SRRid1	SRREXAMPLE1 SRRID2	SRREXAMPLE2"
							    )

							    echo "Starting fastq-dump for all SRR files..."

							    # Step 1: Run fastq-dump for all SRR files sequentially
							    for SRR in "${SRR_FILES[@]}"; do
								        echo "Processing $SRR..."

									    # Check if file exists with or without .sra extension
									        if [[ -f "${SRA_DIR}/${SRR}.sra" ]]; then
											        SRA_FILE="${SRA_DIR}/${SRR}.sra"
												    elif [[ -f "${SRA_DIR}/${SRR}" ]]; then
													            SRA_FILE="${SRA_DIR}/${SRR}"
														        else
																        echo "ERROR: SRA file for ${SRR} not found!" >&2
																	        continue
																		    fi

																		        # Run fastq-dump
																			    fastq-dump --split-files "$SRA_FILE" --outdir "$OUTPUT_DIR"
																		    done

																		    echo "FASTQ dumping complete. Now computing MD5 checksums in pairs..."

																		    # Step 2: Compute MD5 checksums in pairs
																		    rm -f $MD5_FILE  # Remove old checksum file if it exists

																		    for ((i = 0; i < ${#SRR_FILES[@]}; i+=2)); do
																			        SRR1=${SRR_FILES[$i]}
																				    SRR2=${SRR_FILES[$i+1]}

																				        echo "Computing MD5 for $SRR1 and $SRR2..."

																					    # Compute MD5 for both SRRs
																					        md5sum ${OUTPUT_DIR}/${SRR1}_1.fastq >> $MD5_FILE
																						    md5sum ${OUTPUT_DIR}/${SRR1}_2.fastq >> $MD5_FILE
																						        md5sum ${OUTPUT_DIR}/${SRR2}_1.fastq >> $MD5_FILE
																							    md5sum ${OUTPUT_DIR}/${SRR2}_2.fastq >> $MD5_FILE

																							        # Compare checksums
																								    echo "MD5 comparison for $SRR1 and $SRR2:"
																								        diff <(md5sum ${OUTPUT_DIR}/${SRR1}_1.fastq ${OUTPUT_DIR}/${SRR2}_1.fastq) <(md5sum ${OUTPUT_DIR}/${SRR1}_2.fastq ${OUTPUT_DIR}/${SRR2}_2.fastq)

																									    echo "Completed MD5 check for $SRR1 and $SRR2."
																								    done

																								    echo "FASTQ extraction and MD5 pairwise checking completed."

