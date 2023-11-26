#!/bin/bash

partition_file_dir="/home2/yxwu/convert_iqtree_model_to_mrbayes_format/StegnogrammaF_cpDNA_partition.txt"

#處理partition txt變成list
sed -e 's/DNA,\ //g' -e 's/\ =\ .*$//g' *partition.txt > partition_list.txt

rm mb_format_with8space.txt
rm mb_format.txt


#建立字典
declare -A mb_model
mb_model=(
["SYM"]='lset applyto=() nst=6;cj04prset applyto=() statefreqpr=fixed(equal);' 
["SYM+G4"]='lset applyto=() nst=6 rates=gamma;cj04prset applyto=() statefreqpr=fixed(equal);' 
["SYM+I"]='lset applyto=() nst=6 rates=propinv;cj04prset applyto=() statefreqpr=fixed(equal);' 
["SYM+I+G4"]='lset applyto=() nst=6 rates=invgamma;cj04prset applyto=() statefreqpr=fixed(equal);'
["F81+F"]='lset applyto=() nst=1;' 
["F81+F+G4"]='lset applyto=() nst=1 rates=gamma;' 
["F81+F+I"]='lset applyto=() nst=1 rates=propinv;' 
["F81+F+I+G4"]='lset applyto=() nst=1 rates=invgamma;' 
["GTR+F"]='lset applyto=() nst=6;' 
["GTR+F+G4"]='lset applyto=() nst=6;' 
["GTR+F+I"]='lset applyto=() nst=6 rates=propinv;' 
["GTR+F+I+G4"]='lset applyto=() nst=6 rates=invgamma;' 
["HKY+F"]='lset applyto=() nst=2;' 
["HKY+F+G4"]='lset applyto=() nst=2 rates=gamma;' 
["HKY+F+I"]='lset applyto=() nst=2 rates=propinv;' 
["HKY+F+I+G4"]='lset applyto=() nst=2 rates=invgamma;' 
["JC"]='lset applyto=() nst=1;cj04prset applyto=() statefreqpr=fixed(equal);' 
["JC+G4"]='lset applyto=() nst=1 rates=gamma;cj04prset applyto=() statefreqpr=fixed(equal);' 
["JC+I"]='lset applyto=() nst=1 rates=propinv;cj04prset applyto=() statefreqpr=fixed(equal);' 
["JC+I+G4"]='lset applyto=() nst=1 rates=incgamma;cj04prset applyto=() statefreqpr=fixed(equal);' 
["K2P"]='lset applyto=() nst=2;cj04prset applyto=() statefreqpr=fixed(equal);' 
["K2P+G4"]='lset applyto=() nst=2 rates=gamma;cj04prset applyto=() statefreqpr=fixed(equal);' 
["K2P+I"]='lset applyto=() nst=2 rates=propinv;cj04prset applyto=() statefreqpr=fixed(equal);' 
["K2P+I+G4"]='lset applyto=() nst=2 rates=invgamma;cj04prset applyto=() statefreqpr=fixed(equal);'
)


get_line_of_str_infile(){
	#${1} partition_name
	#${2} file_name
	grep -n -w ${1} ${2} | cut -f1 -d:
	return $?
}

#get_partition_name_line
#grep -n rbcL_pos1 StegnogrammaFF2_cpDNA.model | cut -f1 -d:
#這行也一樣awk '/rbcL_pos1/ {print FNR}' StegnogrammaFF2_cpDNA.model

#處理model
index=1
while read sample_name
do
  # echo ${index}
  best_model_BIC_line_tmp=$(grep -n best_model_BIC *.model | cut -f1 -d:)
  best_model_BIC_line=[]
  best_model_BIC_line=(${best_model_BIC_line_tmp// / })
  # echo ${sample_name} 
  partition_name_line=$(get_line_of_str_infile ^${sample_name}: *.model)
  # echo ${partition_name_line}
  
  best_model_BIC_line[${#best_model_BIC_line[@]}]=${partition_name_line}
  
  
  for ((i=1;i<${#best_model_BIC_line[*]};i++))
  do
      for ((a=0;a<${#best_model_BIC_line[*]}-i;a++))
      do
        if [[ ${best_model_BIC_line[$a]} -gt ${best_model_BIC_line[$a+1]} ]]
        then
          temp=${best_model_BIC_line[$a]}
          best_model_BIC_line[$a]=${best_model_BIC_line[$a+1]}
          best_model_BIC_line[$a+1]=$temp
        fi
      done
  done
  # echo ${best_model_BIC_line[@]}
  
  position_found=$(echo ${best_model_BIC_line[@]/${partition_name_line}//} | cut -d/ -f1 | wc -w | tr -d ' ')+1
  model=$(sed -n "${best_model_BIC_line[${position_found}]}p" *.model | cut -f3 -d' ')
  # echo ${sample_name} "is" 
  echo ${mb_model["${model}"]} | sed -e 's/cj04/\n/g' | sed -e "s/()/(${index})/g" >> mb_format.txt
  unset -v model
  unset -v best_model_BIC_line
  
  ((index=index+1))
  
  
done < partition_list.txt

partition_list=$(cat partition_list.txt | sed -e ':a;N;s/\n/, /g;ta')

sed -e 's/DNA,/        charset/g' -e 's/$/;/g' *partition.txt > mb_format_with8space.txt
((index=index-1))
echo "        partition favored =" ${index} ": " ${partition_list}";" >> mb_format_with8space.txt
echo "        set partition = favored;" >> mb_format_with8space.txt
sed -e 's/lset /        lset /g' -e 's/prset /        prset /g' mb_format.txt >> mb_format_with8space.txt

rm mb_format.txt

echo OK!

