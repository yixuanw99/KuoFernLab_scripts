#!/bin/bash

echo "請貼上.model.gz檔與partiion.txt所在的絕對位址 -> "
read modelfilefir

#確認位址格式對不對
test -d $modelfilefir
while [ $? != 0 ]
do
    echo "請輸入位址正確格式 ->"
    read modelfilefir
    test -d $modelfilefir
done

gzip -dr $modelfilefir

#當前位置為這樣，底下有一個partition.txt跟一個.model將依照這兩個檔案進行model的格式轉換
echo "------------------------------------"
echo "以下顯示的資料夾位址為$modelfilefir"
echo -e "底下的所有partition.txt跟.model如下:\n"
ls -1 $modelfilefir | grep "\.model"
ls -1 $modelfilefir | grep "partition.txt"
echo -e "\n確認請輸入y，離開請輸入n\n"
echo -e "\033[31m !!!注意!!! \033[0m 此位置下若有多個partition.txt或.model將會跑出錯誤結果"
echo "[y/n] ->"
read dircheck

#確認y或n到底有沒有打對
while [ $dircheck != "y" ] && [ $dircheck != "n" ]
do
echo -e "----------\n請打y或n"
echo "[y/n] ->"
read dircheck
if [ $dircheck == "y" ] || [ $dircheck == "n" ]
then
    echo "感謝你終於打對了"
    break
else
    echo "再給你一次機會"
fi
done

#y要做啥n要做啥
if [ $dircheck == "y" ]
then
    echo -e "將依照這兩個檔案進行model的格式轉換"
	sleep .3
elif [ $dircheck == "n" ]
then
    echo -e "----------\n請重新執行並確認輸入的位址及裡面檔案是否正確"
    exit 1
else
    echo "啦啦啦請打y或n啦"
fi

partition_file_dir="/home2/yxwu/convert_iqtree_model_to_mrbayes_format/StegnogrammaF_cpDNA_partition.txt"

#處理partition txt變成list
sed -e 's/DNA,\ //g' -e 's/\ =\ .*$//g' $modelfilefir/*partition.txt > $modelfilefir/partition_list.txt

[[ -f "$modelfilefir/mb_format_with8space.txt" ]] && rm $modelfilefir/mb_format_with8space.txt

[[ -f "$modelfilefir/mb_format.txt" ]] && rm $modelfilefir/mb_format.txt

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
  sp='='
  printf ' '
  printf '\b%1s' "$sp"
  sp=${sp#?}${sp%???}
  # echo ${index}
  best_model_list_BIC_line_tmp=$(grep -n best_model_list_BIC *.model | cut -f1 -d:)
  best_model_list_BIC_line=[]
  best_model_list_BIC_line=(${best_model_list_BIC_line_tmp// / })
  # echo ${sample_name} 
  partition_name_line=$(get_line_of_str_infile ^${sample_name}: *.model)
  # echo ${partition_name_line}
  
  best_model_list_BIC_line[${#best_model_list_BIC_line[@]}]=${partition_name_line}
  
  
  for ((i=1;i<${#best_model_list_BIC_line[*]};i++))
  do
      for ((a=0;a<${#best_model_list_BIC_line[*]}-i;a++))
      do
        if [[ ${best_model_list_BIC_line[$a]} -gt ${best_model_list_BIC_line[$a+1]} ]]
        then
          temp=${best_model_list_BIC_line[$a]}
          best_model_list_BIC_line[$a]=${best_model_list_BIC_line[$a+1]}
          best_model_list_BIC_line[$a+1]=$temp
        fi
      done
  done
  # echo ${best_model_list_BIC_line[@]}
  
  position_found=$(($(echo ${best_model_list_BIC_line[@]/${partition_name_line}//} | cut -d/ -f1 | wc -w | tr -d ' ')+1))
  # echo "position_found: "${position_found}
  # echo "best_model_list_BIC_line: " ${best_model_list_BIC_line[@]}
  best_model_list_BIC=$(sed -n "${best_model_list_BIC_line[${position_found}]}p" *.model | cut -f3- -d' ' | sed -e 's/ /\n/g')
  # echo "best_model_list_BIC: "${best_model_list_BIC}
  # printf '%s\n' "$best_model_list_BIC"
  while IFS= read -r best_model
  do
    # echo "----------"${best_model}
    best_model_output=$(echo ${mb_model["${best_model}"]})
    if [[ -z "${best_model_output}" ]]
	then
	  continue
	else
	  echo ${best_model_output} | sed -e 's/cj04/\n/g' | sed -e "s/()/(${index})/g" >> $modelfilefir/mb_format.txt
	  break
	fi
  done < <(printf '%s\n' "$best_model_list_BIC")
  # echo ${best_model_list_BIC}
  #model=$(sed -n "${best_model_list_BIC_line[${position_found}]}p" *.model | cut -f3 -d' ')
  # echo ${sample_name} "is" 
  unset -v best_model_output
  unset -v best_model_list_BIC
  unset -v best_model_list_BIC_line
  
  ((index=index+1))
  
  
done < $modelfilefir/partition_list.txt

partition_list=$(cat $modelfilefir/partition_list.txt | sed -e ':a;N;s/\n/, /g;ta')

sed -e 's/DNA,/        charset/g' -e 's/$/;/g' *partition.txt > $modelfilefir/mb_format_with8space.txt
((index=index-1))
echo "        partition favored =" ${index} ": " ${partition_list}";" >> $modelfilefir/mb_format_with8space.txt
echo "        set partition = favored;" >> $modelfilefir/mb_format_with8space.txt
sed -e 's/lset /        lset /g' -e 's/prset /        prset /g' $modelfilefir/mb_format.txt >> $modelfilefir/mb_format_with8space.txt

rm $modelfilefir/mb_format.txt

echo OK!

