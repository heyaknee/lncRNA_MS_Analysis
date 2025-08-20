#!/bin/bash
#n to run STAR alignment
run_star_alignment() {
	    local SAMPLE=$1
	        local OUT_DIR=$2
		    local FASTQ_DIR=$3  # New input directory parameter

		        echo "Processing sample: $SAMPLE"

			    # Run STAR with proper quoting and corrected file paths
			        STAR --runThreadN 12 \
					         --genomeDir "path_to_your/STAR_Index" \
						          --readFilesIn "${FASTQ_DIR}/${SAMPLE}_1.fastq" "${FASTQ_DIR}/${SAMPLE}_2.fastq" \
							           --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
								            --outSAMtype BAM SortedByCoordinate \
									             --outSAMstrandField intronMotif \
										              --twopassMode Basic \
											               --sjdbGTFfile "path_to_your/gencode.v47.annotation.gtf" \
												                --readFilesCommand cat \
														         --quantMode GeneCounts
				    
				    # Index the BAM file
				        samtools index "${OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

					    echo "STAR alignment completed for ${SAMPLE}, results saved in $OUT_DIR"
				    }

				    # Function to run StringTie
				    run_stringtie() {
					        local SAMPLE=$1
						    local BAM_FILE=$2
						        local OUTPUT_DIR=$3

							    echo "Running StringTie for sample: $SAMPLE"

							        # Run StringTie
								    stringtie "$BAM_FILE" \
									            -p 8 \
										            -G "path_to_your/gencode.v47.annotation.gtf" \
											            -o "${OUTPUT_DIR}/${SAMPLE}.gtf"

								        echo "StringTie completed for ${SAMPLE}, results saved in $OUTPUT_DIR"
								}

								# Set base paths
								BASE_PATH="path to your base directory"  # Adjust this path
								GENOME_DIR="${BASE_PATH}/STAR_Index"
								ANNOTATION_GTF="${BASE_PATH}/gencode.v47.annotation.gtf"

								# Folders containing the fastq files
								FOLDERS=("PPMS" "Non-MS_control" "RRMS" "SPMS") #these folders should contain the fastq files of different stages of the disease 

								# Create log file for the pipeline
								LOG_FILE="${BASE_PATH}/pipeline_log.txt"
								echo "Pipeline started at $(date)" > "$LOG_FILE"

								# Loop through each folder 
								for folder in "${FOLDERS[@]}"; do
									    # Set the output folder based on the group
									        if [[ "${folder}" == "PPMS" ]]; then
											        OUT_DIR="${BASE_PATH}/${folder}/star_results_ppms"
												        STRINGTIE_OUT_DIR="${BASE_PATH}/${folder}/stringtie_results_ppms"
													    elif [[ "${folder}" == "Non-MS_control" ]]; then
														            OUT_DIR="${BASE_PATH}/${folder}/star_results_nonms"
															            STRINGTIE_OUT_DIR="${BASE_PATH}/${folder}/stringtie_results_nonms"
																        elif [[ "${folder}" == "RRMS" ]]; then
																		        OUT_DIR="${BASE_PATH}/${folder}/star_results_rrms"
																			        STRINGTIE_OUT_DIR="${BASE_PATH}/${folder}/stringtie_results_rrms"
																				    elif [[ "${folder}" == "SPMS" ]]; then
																					            OUT_DIR="${BASE_PATH}/${folder}/star_results_spms"
																						            STRINGTIE_OUT_DIR="${BASE_PATH}/${folder}/stringtie_results_spms"
																							        fi

																								    # Create output directories with proper space handling
																								        mkdir -p "$OUT_DIR"
																									    mkdir -p "$STRINGTIE_OUT_DIR"

																									        # Input directory for this folder
																										    INPUT_DIR="${BASE_PATH}/${folder}"

																										        # Get total number of samples
																											    TOTAL_SAMPLES=$(ls "${INPUT_DIR}/"*_1.fastq | wc -l)

																											        # Process samples with proper quoting and file extension
																												    for ((i=1; i<=TOTAL_SAMPLES; i++)); do
																													            file_1=$(ls "${INPUT_DIR}/"*_1.fastq | head -n $i | tail -n 1)
																														            SAMPLE=$(basename "$file_1" _1.fastq)
																															            
																															            # Verify paired-end file exists
																																            if [ -f "${INPUT_DIR}/${SAMPLE}_2.fastq" ]; then
																																		                echo "Processing sample $i of $TOTAL_SAMPLES: $SAMPLE"
																																				            echo "Processing sample $i of $TOTAL_SAMPLES: $SAMPLE" >> "$LOG_FILE"
																																					                
																																					                # Run STAR alignment
																																							            run_star_alignment "$SAMPLE" "$OUT_DIR" "$INPUT_DIR"
																																								                echo "STAR alignment completed for $SAMPLE" >> "$LOG_FILE"
																																										            
																																										            # Run StringTie
																																											                BAM_FILE="${OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
																																													            run_stringtie "$SAMPLE" "$BAM_FILE" "$STRINGTIE_OUT_DIR"
																																														                echo "StringTie completed for $SAMPLE" >> "$LOG_FILE"
																																																        else
																																																		            echo "ERROR: Missing paired-end file for ${SAMPLE}"
																																																			                echo "ERROR: Missing paired-end file for ${SAMPLE}" >> "$LOG_FILE"
																																																					        fi
																																																						    done
																																																					    done

																																																					    echo "✅ StringTie Transcript Assembly Completed!"
																																																					    echo "Pipeline completed at $(date)" >> "$LOG_FILE"

