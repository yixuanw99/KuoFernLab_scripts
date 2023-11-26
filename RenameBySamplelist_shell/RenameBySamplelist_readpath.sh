#! /bin/bash

echo "請貼上fasta與samplelsit.txt所在的絕對位址 -> "
read fastafiledir

#確認位址格式對不對
test -d $fastafiledir
while [ $? != 0 ]
do
    echo "請輸入位址正確格式 ->"
    read fastafiledir
    test -d $fastafiledir
done

#列出資料夾內容並等待確認
echo "------------------------------------"
echo "以下顯示的資料夾位址為$fastafiledir"
echo -e "請確認此資料夾內容是否正確\n正確請輸入y，錯誤請輸入n"
ls -1 $fastafiledir
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
    echo "將確認fasta是否備齊"
elif [ $dircheck == "n" ]
then
    echo -e "----------\n請重新執行並確認輸入的位址是否正確"
    exit 1
else
    echo "啦啦啦請打y或n啦"
fi

#declare some variable from list
SampleList=$fastafiledir'/SampleList.txt'
sed -i "s/\r//g" $SampleList
Column_Number=($(awk 'NR==1{print NF}' $SampleList))
if [ $? -eq 2 ]
then
    echo -e "------------------------------------\n請確認你的sample list是否命名為此檔名\"SampleList.txt\"，檔名不一樣無法執行喔"
    exit 1
fi
Locus_Number=$((${Column_Number[0]}-1))
Row_Number=($(awk 'END{print FNR}' $SampleList))

#read sample list into matrix
declare -A Sample_list
num_columns=$Column_Number
num_rows=$Row_Number
for ((j=1;j<=num_columns;j++)) do
    for ((i=1;i<=num_rows;i++)) do
        tmp=($(awk -v i=$i -v j=$j 'FNR == i {printf "%s", $j}' $SampleList))
        Sample_list[$j,$i]=$tmp
    done
    unset tmp
done

#確認fasta是否備齊
for ((i=2;i<=num_columns;i++)) do
    echo "searching" ${Sample_list[$i,1]}
    ls $fastafiledir | grep -q "^${Sample_list[$i,1]}\.fasta"
    if [ $? -eq 1 ]
    then
        echo "----->not find ${Sample_list[$i,1]}.fasta--QAQ"
        echo "請補上再重跑一次喔"
        exit 1
    else
        echo "----->find ${Sample_list[$i,1]}.fasta--^o^"
    fi
done

#到這步還沒exit 1的話就開始跑囉
echo "將開始執行重命名腳本"

#replace all space with underline, rename each fasta header by sample list
for ((k=2;k<=num_columns;k++)) do
    InputFasta=$fastafiledir'/'${Sample_list[$k,1]}'.fasta'
	OutputFasta=$fastafiledir'/'${Sample_list[$k,1]}'_renamed.fas'
	FastaUnderLine=$fastafiledir'/'${Sample_list[$k,1]}'_Underline.fasta'
	sed 's/[[:space:]]/_/g' $InputFasta > $FastaUnderLine
	cp $FastaUnderLine $OutputFasta
	for ((l=2;l<=num_rows;l++)) do
        grep -n ${Sample_list[$k,$l]} ${FastaUnderLine}
		if [ $? -ne 0 ] ;then
			echo "not find accession number(${Sample_list[$k,$i]}) in input fasta"
		else
			sed -i '/'${Sample_list[$k,$l]}'/c \>'${Sample_list[1,$l]}'' ${OutputFasta}
			echo "processing locus ${Sample_list[$k,1]}"
			echo "processing sample ${Sample_list[1,$l]}"
			echo "replace success: find accession number(${Sample_list[$k,$l]}) in input fasta file and replaced by(>${Sample_list[1,$l]})"
			echo "------------------------------------------------------------"
		fi
    done
    rm $FastaUnderLine
done

#sample list要把名稱空格全部取代成底線
#fasta也要把名稱空格全部取代成底線(已做好)
#fasta要存成單一換行格式的(bioedit的不行，aliview的才可以)

# 將sample list讀入並建立待處理的matrix
# 用sample list檢查fasta是否備齊(之後做)
# 從第一個locus開始做(sample list的第二列)
# 	sample list的每一行(每個sample)重複做
# 		在陣列中第一列(header)中搜尋sample list的accession number
#		找不到的話就說找不到(建立missing sequence_之後做)
# 		找到的話就重命名那行為sample list的taxa
# 	一路搜尋到最後一個accession number
# 做到最後一個locus(sample list的最後一列)