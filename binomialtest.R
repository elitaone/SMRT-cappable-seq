# Created by Bo Yan, released in 2017
# used for performing the binomial test to determine the significant 3'end
# usage: R binomialtest.R --no-save input.bed output.bed 0.2 < binomialtest.R
# See explaination in README.txt

args = commandArgs()
input = args[4]
output = args[5]
percent = as.numeric(args[6])
print(input)
print(output)
print(percent)

new = read.csv(input, sep='\t', header=F)
times=1000*percent
## since I use >=10 reads as cutoff for 3'end determination, I can remove the entries that have less than 10 for the same reads to reduce the amount of caculation
new = new[which(new$V8>=10),] 
## caculate the distance, which is the length of the largest transcript for all the transripts with the same TSS
positive = new[which(new$V6=='+'),]
positive$distance = positive$V4-positive$V2
negative = new[which(new$V6=='-'),]
negative$distance = negative$V3-negative$V4
new = rbind(positive, negative)

## perform binomial test for each entry and add the pvalue to new$pvalue, and add tag for 95% and 99% confidence
new$pvalue = 1
new$oneside95 = 0
new$oneside99 = 0
for(i in 1:nrow(new)){
  prob = min(1/new[i,9]*times, percent)
  new[i,10] = 1- pbinom(new[i,8],new[i,7],prob)
  if (new[i,10]<=0.05) {
    new[i,11]=1
  } 
  if (new[i,10]<=0.01) {
    new[i,12]=1
  }
}

new=new[which(new$oneside95==1),]
write.table(new, output, sep='\t', quote=F, row.names=F, col.names=F)
q()