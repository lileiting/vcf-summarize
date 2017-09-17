#!/bin/bash

args()
{	echo "Summarize and plot(auto detect variable type) for all or specified variables in a vcf file (except CHROM POS ID fields,which I think not necessary).
Extract arbitrary fixed fields or values in sample field to TAB delimited file. 
Natural support for multi-sample vcf

before using
chmod +x vcf-summarize.sh

Simplest usage: vcf-summarize.sh -f filename.vcf -a #extract and summarize all fields and subfields 
For large file: nohup vcf-summarize.sh -f filename.vcf -a &

Usage 
# -f [required] Take 1 file. The target vcf file. Support plain txt and gz,bz file
# -a [optional] Take 0 argument. If specified, extract and summarize all variables
# -q [optional] Take 0 argument. If specified, will skip REF, ALT, QUAL, FILTER fields
# -c [optional] Take 1 string. e.g. -c \"chr1\" will limit analysis to records with CHROM fields equal to chr1 
# -i [optional] Take 1 string. e.g. -i \"CHROM POS\" will used CHR_POS [default] as the index of extracted tables. You can choose any combination of \"CHROM POS ID REF ALT\". The idea is to generate unique index with the smallest number of fields.   
# -I [optional] Take 1 string. e.g. -I \"AN DB\" will extract and summarize AN and DB subfields in INFO field. Will overwrite option -a, which analyze all subfields in INFO.
# -F [optional] Take 1 string. e.g. -F \"GT AD DP\" will extract and summarize GT AD DP subfields in sample columns. Will overwrite option -a, which analyze all subfields in sample columns
# -s [optional] Take 0 argument. If specified, just do data extraction. Suppress summarization and plotting
# -o [optional] Take 1 string. The output directory name. Default is vcfsummarize
# -h [optional] show this help

#This is a personal effort without any funding support
#Report suggestion and bug to ruansun@163.com";
	exit;
}

log1="vcfsummary.log";
echo see complete log at $log1;

allfields="";
emitRAQF=1;
chrome="";
index="CHROM POS";
inputfile="";
infovar="";
formvar="";
suppressplot=0;
outputdir="vcfsummarize";
while getopts f:aqc:i:I:F:so:h opt
do
	case $opt in
	f) inputfile="$OPTARG";;
	a) allfields=1;;
	q) emitRAQF=0;;
	c) chrome="$OPTARG";;
	i) index="$OPTARG";;	
	I) infovar="$OPTARG";;
	F) formvar="$OPTARG";;
	s) suppressplot=1;;
	o) outputdir="$OPTARG";;
	h) args;;
	esac
done

if [[ -d "$outputdir" ]]
then
	echo "$outputdir exists, wish to overwrite?[Y/N]";
	read answer;
	if [[ "$answer" == Y* || "$answer" == y* ]]
	then
		echo "overwrite $outputdir";
		mkdir -p "$outputdir";
	else
		exit;
	fi
else
	mkdir -p "$outputdir";
fi 

echo write output to $outputdir;

if [[ -z $inputfile ]]
then
	echo '-f [inputfile] missing';
	args;
fi
 
filename=${inputfile##*/};

#test if input is ZIP file
zipfile=`file $inputfile | grep 'zip' | wc -l`; 
if [[ $zipfile -eq 0 ]]
then
	streamer="cat";
else
	streamer="gunzip -c";
	filename=${filename%.*};
fi



eval `$streamer $inputfile | awk -F'\t' '/^#[^#]/ {for(i=10;i<=NF;i++){smpname=smpname":"$i;}exit;} 
	END{print "smpname="smpname;}'` ;

smpname=${smpname#:};
smpnum=`echo $smpname | tr ':' '\n' | wc -l`;
echo found $smpnum samples | tee -a  $log1;

echo parsing variable type | tee -a  $log1; 
#find non numeric variable
eval `$streamer $inputfile | head -2000 | awk -F"\t" '/^##INFO|^##FORMAT/ {var=$1;
	sub(/,.*/,"",var);
	sub(/.*=/,"",var);
	if($1 ~ /^##INFO/)
	{	infovar=infovar"-"var;
	}
	if($1 ~ /^##FORMAT/)
	{	formvar=formvar"-"var;
	}
	if($1 !~ /Number=1,/ || $1 ~ /Type=String/)
	{	categ=categ"-"var;
	}
	}
	END{print "categ="categ;
	print "allinfovar="infovar;
	print "allformvar="formvar;}'`;

tmpcateg="`echo $categ | tr '-' ' '`"; 
echo auto detect REF ALT FILTER $tmpcateg as categorical variable | tee -a  $log1; 

if [[ -z $infovar && $allfields -eq 1 ]]
then
	allinfovar="${allinfovar#-}";
	infovar="`echo $allinfovar | tr '-' ' '`";
fi

if [[ -z $formvar && $allfields -eq 1 ]]
then
	allformvar="${allformvar#-}";
	formvar="`echo $allformvar | tr '-' ' '`";
fi

if [[ ! -z $chrome ]]
then 
	echo process chromosome $chrome | tee -a $log1;
else
	echo process all chromosomes in $inputfile | tee -a $log1;
fi

if [[ $emitRAQF -eq 1 ]]
then
	echo process REF ALT QUAL FILTER fields | tee -a $log1;
else
	echo skip REF ALT QUAL FILTER fields | tee -a $log1;
fi

echo process INFO fields \"$infovar\" | tee -a  $log1;
echo process sample fields \"$formvar\" | tee -a  $log1;
echo output fixed field extracts to $outputdir"/""$filename".extract | tee -a  $log1;
echo output sample field extracts to $outputdir"/""$filename".VARNAME | tee -a $log1;
echo output sample field summary to $outputdir"/"outputforsamples | tee -a  $log1; 
 
#generate data extracts
$streamer $inputfile | awk -F'\t' -v cusindex="$index" -v chrome="$chrome" -v outputdir="$outputdir" -v smpname="$smpname" -v filename=$filename -v emitRAQF=$emitRAQF  -v infovar="$infovar" -v formvar="$formvar" 'BEGIN{OFS="\t";
	cusindexhash["CHROM"]=1;
	cusindexhash["POS"]=2;
	cusindexhash["ID"]=3;
	cusindexhash["REF"]=4;
	cusindexhash["ALT"]=5;
	cusindexarrsize=split(cusindex,cusindexarr," ");
	iter=1;
	for(i=1;i<=cusindexarrsize;i++)
	{	if(cusindexarr[i] in cusindexhash)
		{	cusindexid[iter]=cusindexhash[cusindexarr[i]];
			iter+=1;
		}
	}
	if(length(cusindexid)==0)
	{	cusindexid[1]=1;
		cusindexid[2]=2;
	}
	if(emitRAQF==1)
	{	printf "REF\tALT\tQUAL\tFILTER\t";
	}
	infoarrsize=split(infovar,infoarr," ");
	for(i=1;i<=infoarrsize;i++)
	{	printf "%s\t",infoarr[i];
	}
	printf "\n";
	formarrsize=split(formvar,formarr," ");
	gsub(/:/,"\t",smpname);
	for(i=1;i<=formarrsize;i++)
	{	printf "INDEX\t%s\n",smpname > outputdir"/"filename"."formarr[i];
	}
	}
	
	 /^[^#]/ {
	if((chrome != "") && ($1 != chrome))
	{	next;
	}
	if($7 == ".")
	{	$7="Not_Set";
	}	
	if(emitRAQF==1)
	{	printf "%s\t%s\t%s\t%s\t",$4,$5,$6,$7;
	}
	size=split($8,arr,";");
	for(i=1;i<=size;i++)
	{	split(arr[i],subarr,"=");
		if(2 in subarr)
		{	rec[subarr[1]]=subarr[2];
		}
		else
		{	rec[subarr[1]]=subarr[1];
		}
	}
	for(i=1;i<=infoarrsize;i++)
	{	if(infoarr[i] in rec)
		{	printf "%s\t",rec[infoarr[i]];
		}
		else
		{	printf "%s\t","NA";
		}
	}
	printf "\n";
	size=split($9,arr,":");
	for(i=1;i<=size;i++)
	{	formhash[arr[i]]=i;
	}

	for(j=1;j<=formarrsize;j++)
	{	sep="_";
		for(idx=1;idx<=length(cusindexid);idx++)
		{	if(idx==length(cusindexid))
			{	sep="\t";
			}
			printf "%s%s",$cusindexid[idx],sep > outputdir"/"filename"."formarr[j];
		}
	}
	for(i=10;i<=NF;i++)
	{	size=split($i,smp,":");
		sep="\t";
		if(i==NF)
		{	sep="";
		}
		for(j=1;j<=formarrsize;j++)
		{	
			if(formarr[j] in formhash && formhash[formarr[j]] in smp)
			{	printf "%s%s",smp[formhash[formarr[j]]],sep > outputdir"/"filename"."formarr[j];
			}
			else
			{	printf "NA%s",sep > outputdir"/"filename"."formarr[j];
			}
		}
	}
	split("",formhash);	#clear content in formhash. this is more portable than delete formhash
	for(j=1;j<=formarrsize;j++)
	{	printf "\n" > outputdir"/"filename"."formarr[j];
	}
	}' > "$outputdir"/"$filename".extract;

if [[ $suppressplot -eq 0 ]]
then
linecheck=`wc -l "$outputdir"/"$filename".extract | cut -d ' ' -f 1`; 
if [[ $emitRAQF -eq 1 && $linecheck -le 1 ]]
then
	echo no record found;
	exit;
fi

	if [[ $linecheck -gt 1 ]]
	then
Rscript <(echo 'args=commandArgs(TRUE);
	filename=args[1];
	catvar=strsplit(args[2],"_")[[1]];
	catvar=c(catvar,"REF","ALT","FILTER");
	outputdir=args[3];
	print("start R session for summarizing and plotting fixed field variable");
	inputfile=paste(outputdir,"/",filename,".extract",sep="");
	cat("read input from ",inputfile,"\n");
	tmp=read.table(inputfile,header=T,sep="\t");
	tmp=tmp[,1:(ncol(tmp)-1)];
	tr=lapply(colnames(tmp),function(xi){
		if(! xi %in% catvar)
		{	cat("convert ",xi," to numeric\n");
			tmp[,xi]<<-as.numeric(as.character(tmp[,xi]));
		}});
	tr=lapply(colnames(tmp),function(xi){
		outputfile=paste(outputdir,"/summary.",filename,sep="");
		cat("write summary for fixed fields ",xi," to ",outputfile,"\n");
		res=summary(tmp[,xi],maxsum=50);
		write.table(paste(">>",xi,sep=""),file=outputfile,quote=F,col.names=F,row.names=F,sep="\t",append=T);
		write.table(as.matrix(res),file=outputfile,quote=F,col.names=F,sep="\t",append=T);
		cat("generate plot for ",xi,"\n");
		png(paste(outputdir,"/",filename,".",xi,".png",sep=""),height=1500,width=750);
		par(mfrow=c(2,1),mar=c(4,6,4,4));
		if(! xi %in% catvar)
		{	hist(tmp[,xi],main=xi,xlab="",cex.axis=2,cex.lab=2,col="grey");
			boxplot(tmp[,xi],main=xi,cex.axis=2,cex.lab=2,col="grey");
		}
		else
		{	res=summary(tmp[,xi],maxsum=20);
			pie(res,main=xi,cex.main=2,cex=2);
			barplot(res,col=c(1:length(res)),main=xi,xlab="",cex.main=2,cex.axis=2);
		}
		dev.off();
	})') "$filename" "$categ" "$outputdir" >> $log1;
	fi

	if [[ ! -z $formvar ]]
	then
mkdir -p "$outputdir"/outputforsamples;
Rscript <(echo 'args=commandArgs(TRUE);
	filename=args[1];
	formvar=strsplit(args[2]," ")[[1]];
	catvar=strsplit(args[3],"-")[[1]];
	outputdir=args[4];
	print("start R session for summarizing and plotting sample field variable");
	tr=lapply(formvar,function(xi){
		inputfile=paste(outputdir,"/",filename,".",xi,sep="");
		cat("read input from ",inputfile,"\n");
		tmp=read.table(inputfile,header=T,sep="\t");
		tmp=tmp[,2:ncol(tmp)];
		if(xi %in% catvar)
		{	tr=lapply(1:ncol(tmp),function(j){tmp[,j]<<-as.factor(tmp[,j]);});
		}
		else
		{	tr=lapply(1:ncol(tmp),function(j){tmp[,j]<<-as.numeric(as.character(tmp[,j]));});
		}
		tr=lapply(colnames(tmp),function(yi){
			cat("generate summary for ",xi," attributes of sample ",yi,"\n");
			res=summary(tmp[,yi],maxsum=50);
			write.table(paste(">>",yi,sep=""),file=paste(outputdir,"/outputforsamples/summary.",filename,".",xi,sep=""),sep="\t",quote=F,col.names=F,row.names=F,append=T);
			write.table(as.matrix(res),file=paste(outputdir,"/outputforsamples/summary.",filename,".",xi,sep=""),sep="\t",quote=F,col.names=F,append=T);
			cat("generate plot for ",xi," attributes of sample ",yi,"\n");
			png(paste(outputdir,"/outputforsamples/",filename,".",xi,".",yi,".png",sep=""),height=1500,width=750);
			par(mfrow=c(2,1),mar=c(4,6,4,4));
			if(! xi %in% catvar)
			{	if(sum(!is.na(tmp[,yi]))>0)
				{	hist(tmp[,yi],main=paste(yi,xi,sep=" "),xlab="",cex.axis=2,cex.lab=2,col="grey");
					boxplot(tmp[,yi],main=paste(yi,xi,sep=" "),cex.axis=2,cex.lab=2,col="grey");
				}
				else
				{	cat("WARNING:No value exists for ",xi," attribute of sample ",yi,"\n");
				}
			}
			else
			{	res=summary(tmp[,yi],maxsum=20);
				pie(res,main=paste(yi,xi,sep=" "),cex.main=2,cex=2);
				barplot(res,col=c(1:length(res)),main=paste(yi,xi,sep=" "),xlab="",cex=2,cex.main=2,cex.axis=2);
			};
			dev.off();
			})
	})') "$filename" "$formvar" "$categ" "$outputdir" >> $log1 ;


fi
fi



