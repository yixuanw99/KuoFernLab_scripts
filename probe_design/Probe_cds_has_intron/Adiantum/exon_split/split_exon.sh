#! /bin/bash


#remove previous temporary file
[[ -f "exon_split.fas" ]] && rm exon_split.fas

#input:KAI_cds_from_genomic.fna_oneline.fas all cds sequence
grep -v "complement" --no-group-separator KAI_cds_oneline.fas | grep -A1 "join" --no-group-separator | sed 's/.*\[protein_id=/>/;s/\] \[location.*join/ join/g' | sed -r 's/join\(//g;s/\)+\].*//g' > containing_intron_sequence.fas
awk '/KAI/{print $2}' containing_intron_sequence.fas > location_of_line_containKAI.txt
awk -F, '{for(i=1; i<= NF; i++){print $i > "location_for"NR".txt"};printf("\n");}' location_of_line_containKAI.txt
sed -i 's/\.\./ /g;s/>//g;s/<//g' location_for*.txt

header_of_line_containKAI=($(awk '/KAI/{print $1}' containing_intron_sequence.fas))
echo "header_of_line_containKAI is "${header_of_line_containKAI[@]}
intron_sequence_amount=$(awk 'BEGIN{ i=0 } /KAI/{i++} END{print i}' containing_intron_sequence.fas)

for ((i=1; i<=$intron_sequence_amount; i++))
do
	exon_start=($(awk '{print $1}' location_for${i}.txt))
	exon_end=($(awk '{print $2}' location_for${i}.txt))
	echo "-----sequence "$i" calculating exon length-----"
	[[ ${#exon_start[@]} -ne ${#exon_end[@]} ]] && exon_end[${#exon_start[@]}-1]=${exon_start[-1]}
	echo "exon start is "${exon_start[@]}
	echo "exon end is "${exon_end[@]}
	for ((j=0; j<${#exon_start[@]}; j++))
	do
		echo "exon_end "$j" is "${exon_end[j]}
		echo "exon_start "$j" is "${exon_start[j]}
		exon_length+=("$((${exon_end[j]}-${exon_start[j]}+1))")
		echo "===>exon length is "${exon_length[j]}
	done
	echo "#####sequence "$i" calculating intron length#####"
	for ((k=0; k<${#exon_start[@]}-1; k++))
	do
		echo "exon_start "$k+1" is "${exon_start[k+1]}
		echo "exon_end "$k" is "${exon_end[k]}
		intron_length+=("$((${exon_start[k+1]}-${exon_end[k]}-1))")
		echo "===>intron length is "${intron_length[k]}
	done
	# grep -A1 "${header_of_line_containKAI[i-1]}" containing_intron_sequence.fas >> exon_split_for${i}.fas
	for ((l=0; l<${#exon_start[@]}; l++))
	do
		last_element=$((${#exon_start[@]}-1))
		if [[ $l -eq 0 ]];then
		  printf "%s_exon%d_intronlen_left0_right%s\n" "${header_of_line_containKAI[i-1]}" "$((${l}+1))" "${intron_length[l]}" >> exon_split_for${i}.fas
		elif [[ $l -eq $last_element ]];then
		  printf "%s_exon%d_intronlen_left%s_right0\n" "${header_of_line_containKAI[i-1]}" "$((${l}+1))" "${intron_length[l-1]}" >> exon_split_for${i}.fas
		else
		  printf "%s_exon%d_intronlen_left%s_right%s\n" "${header_of_line_containKAI[i-1]}" "$((${l}+1))" "${intron_length[l-1]}" "${intron_length[l]}" >> exon_split_for${i}.fas
		fi
		# printf "%s_exon%d_intronlen_l%sr%\n" "${header_of_line_containKAI[i-1]}" "$((${l}+1))" "${intron_length[l]}" >> exon_split_for${i}.fas
		if [ $l -eq 0 ];then
			position1="1"
			position2=${exon_length[l]}
		else
			position1=$((${position2}+1))
			position2=$((${position2}+${exon_length[l]}))
		fi
		grep -A1 "${header_of_line_containKAI[i-1]}" containing_intron_sequence.fas | tail -n +2 | cut -c${position1}-${position2}  >> exon_split_for${i}.fas
	done
	unset exon_length
	unset intron_length
done

cat exon_split_for*.fas > exon_split.fas
rm exon_split_for*.fas
rm location_for*.txt