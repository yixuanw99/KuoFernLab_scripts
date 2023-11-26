#! /bin/bash

echo "program is running: add intron by \"n\""

#remove previous temporary file
[[ -f "intron_added.fas" ]] && rm intron_added.fas

#input:Adiantum_capillusveneris_annotation_oneline.fas all cds sequence
grep -v "complement" --no-group-separator Adiantum_capillusveneris_annotation_oneline.fas | grep -A1 "join" --no-group-separator | sed 's/.*\[protein_id=/>/;s/\] \[location.*join/ join/g' | sed -r 's/join\(//g;s/\)+\].*//g' > containing_intron_sequence.fas
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
	grep -A1 "${header_of_line_containKAI[i-1]}" containing_intron_sequence.fas >> intron_added_for${i}.fas
	for ((l=0; l<${#exon_start[@]}-1; l++))
	do
		#here
		gapnumber=${intron_length[l]}
		if [ $l -eq 0 ];then
			posision=${exon_length[l]}
		else
			posision=$((${posision}+${exon_length[l]}+${intron_length[l-1]}))
		fi
		echo "adding intron: "$gapnumber"gaps added into sequence posision: "$posision
		sed -i "/${header_of_line_containKAI[i-1]}/{n;s/./&$(printf '%0.s-' $(seq 1 $gapnumber))/$posision}" intron_added_for${i}.fas
	done
	unset exon_length
	unset intron_length
done

cat intron_added_for*.fas > intron_added.fas
rm intron_added_for*.fas
rm location_for*.txt