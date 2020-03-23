#Tables should contain DNA sequences in the first columns, and mapped taxonomy path in the following columns.
#Sort Tables by Alphabetical Order of the Sequences
table1 <- table1[order(table1$Sequences),]
#Substitute NA with "NA" string
table1$row1<-as.character(table1$row1)
table1$row2<-as.character(table1$row2)
table1$row3<-as.character(table1$row3)
table1$row4<-as.character(table1$row4)
table1$row5<-as.character(table1$row5)
table1$row6<-as.character(table1$row6)
table1$row7<-as.character(table1$row7)
table1$row8<-as.character(table1$row8)
table1$row1[is.na(table1$row1)]<-"NA"
table1$row2[is.na(table1$row2)]<-"NA"
table1$row3[is.na(table1$row3)]<-"NA"
table1$row4[is.na(table1$row4)]<-"NA"
table1$row5[is.na(table1$row5)]<-"NA"
table1$row6[is.na(table1$row6)]<-"NA"
table1$tow7[is.na(table1$tow7)]<-"NA"
table1$tow8[is.na(table1$row8)]<-"NA"

#unify table column names
row.names(table1)<-NULL
colnames(table1)<-c("Sequence","pr2_r1","pr2_r2","pr2_r3","pr2_r4","pr2_r5","pr2_r6","pr2_r7","pr2_r8")

#Two Table Comparison: Compare two tables and record how many sequences have exactly the same mapping results, how many have different names, how many have one NA and one name and vice versa.
#All Agree
Agree<-c()
for(i in 1:nrow(table1){
  if(table1[i,1]==table2[i,1] && table1[i,2]==table2[i,2] && table1[i,3]==table2[i,3] && table1[i,4]==table2[i,4] && table1[i,5]==table2[i,5] && table1[i,6]==table2[i,6] &&table1[i,7]==table2[i,7] &&LCAS2P[i,8]==LCAP2P[i,8] &&LCAS2P[i,9]==LCAP2P[i,9]){
    Agree<-append(Agree,i)
  }
}

length(Agree)
#Different Names
Name_Name<-c()
for(i in 1:nrow(table1){
  for(j in 1:ncol(table1){
    if(table1[i,j]!=table2[i,j] && table1[i,j]!="NA" && table2[i,j]!="NA"){
      Name_Name<-append(Name_Name,i)
    }
  }
}
length(unique(Name_Name))

#table1 NA, table 2 name
NA_Tax<-c()
for(i in 1:nrow(table1){
  for(j in 1:ncol(table1){
    if(table1[i,j]=="NA" && table2[i,j]!="NA"){
      NA_Tax<-append(NA_Tax,i)
    }
  }
}
length(unique(NA_Tax))

#table1 name, table2 NA
Tax_NA<-c()
for(i in 1:nrow(table1){
  for(j in 1:ncol(table1){
    if(table1[i,j]!="NA" && table2[i,j]=="NA"){
      Tax_NA<-append(Tax_NA,i)
    }
  }
}
length(unique(Tax_NA))

#table1 NA, table2 NA
NA_NA<-c()
for(i in 1:nrow(table1){
  for(j in 1:ncol(table1){
    if(table1[i,j]=="NA" && table2[i,j]=="NA"){
      NA_NA<-append(NA_NA,i)
      break
    }
  }
}
length(unique(NA_NA))

 #create a barplot based on the results above
      
bars1<-c(length(unique(Agree)),length(unique(NA_Tax)),length(unique(Tax_NA)),length(unique(Name_Name)))
lbls1<-c("All Agree","NAinTable1","NAinTable2","Different Name")
barplot(bars1, xlab="Conditions",names.arg=lbls1, ylab = "Number of Occurance",ylim=c(0,20000),col=rainbow(length(bars1)), main="Table1 vs. Table2")

      
#barplots for taxonomy ranks: A grouped barplot is created here to illustrate the number of mapped results for each taxonomy rank in each table.
Rank_Comparison<-function(table1,table2,table3,table4,table5,table6,table7,table8){
  t1r1=0
  t2r1=0
  t3r1=0
  t4r1=0
  t5r1=0
  t6r1=0
  for(i in 1:nrow(table1)){
    if(table1[i,2]!="NA"){
      t1r1=t1r1+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,2]!="NA"){
      t2r1=t2r1+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,2]!="NA"){
      t3r1=t3r1+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,2]!="NA"){
      t4r1=t4r1+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,2]!="NA"){
      t5r1=t5r1+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,2]!="NA"){
      t6r1=t6r1+1
    }
  }
  t1r2=0
  t2r2=0
  t3r2=0
  t4r2=0
  t5r2=0
  t6r2=0
  for(i in 1:nrow(table1)){
    if(table1[i,3]!="NA"){
      t1r2=t1r2+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,3]!="NA"){
      t2r2=t2r2+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,3]!="NA"){
      t3r2=t3r2+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,3]!="NA"){
      t4r2=t4r2+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,3]!="NA"){
      t5r2=t5r2+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,3]!="NA"){
      t6r2=t6r2+1
    }
  }
  t1r3=0
  t2r3=0
  t3r3=0
  t4r3=0
  t5r3=0
  t6r3=0
  for(i in 1:nrow(table1)){
    if(table1[i,4]!="NA"){
      t1r3=t1r3+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,4]!="NA"){
      t2r3=t2r3+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,4]!="NA"){
      t3r3=t3r3+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,4]!="NA"){
      t4r3=t4r3+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,4]!="NA"){
      t5r3=t5r3+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,4]!="NA"){
      t6r3=t6r3+1
    }
  }
  t1r4=0
  t2r4=0
  t3r4=0
  t4r4=0
  t5r4=0
  t6r4=0
  for(i in 1:nrow(table1)){
    if(table1[i,5]!="NA"){
      t1r4=t1r4+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,5]!="NA"){
      t2r4=t2r4+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,5]!="NA"){
      t3r4=t3r4+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,5]!="NA"){
      t4r4=t4r4+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,5]!="NA"){
      t5r4=t5r4+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,5]!="NA"){
      t6r4=t6r4+1
    }
  }
  t1r5=0
  t2r5=0
  t3r5=0
  t4r5=0
  t5r5=0
  t6r5=0
  for(i in 1:nrow(table1)){
    if(table1[i,6]!="NA"){
      t1r5=t1r5+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,6]!="NA"){
      t2r5=t2r5+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,6]!="NA"){
      t3r5=t3r5+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,6]!="NA"){
      t4r5=t4r5+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,6]!="NA"){
      t5r5=t5r5+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,6]!="NA"){
      t6r5=t6r5+1
    }
  }
  t1r6=0
  t2r6=0
  t3r6=0
  t4r6=0
  t5r6=0
  t6r6=0
  for(i in 1:nrow(table1)){
    if(table1[i,7]!="NA"){
      t1r6=t1r6+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,7]!="NA"){
      t2r6=t2r6+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,7]!="NA"){
      t3r6=t3r6+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,7]!="NA"){
      t4r6=t4r6+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,7]!="NA"){
      t5r6=t5r6+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,7]!="NA"){
      t6r6=t6r6+1
    }
  }
  t1r7=0
  t2r7=0
  t3r7=0
  t4r7=0
  t5r7=0
  t6r7=0
  for(i in 1:nrow(table1)){
    if(table1[i,8]!="NA"){
      t1r7=t1r7+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,8]!="NA"){
      t2r7=t2r7+1
    }
  }
  for(i in 1:24126){
    if(table3[i,8]!="NA"){
      t3r7=t3r7+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,8]!="NA"){
      t4r7=t4r7+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,8]!="NA"){
      t5r7=t5r7+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,8]!="NA"){
      t6r7=t6r7+1
    }
  }
  t1r8=0
  t2r8=0
  t3r8=0
  t4r8=0
  t5r8=0
  t6r8=0
  for(i in 1:nrow(table1)){
    if(table1[i,9]!="NA"){
      t1r8=t1r8+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table2[i,9]!="NA"){
      t2r8=t2r8+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table3[i,9]!="NA"){
      t3r8=t3r8+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table4[i,9]!="NA"){
      t4r8=t4r8+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table5[i,9]!="NA"){
      t5r8=t5r8+1
    }
  }
  for(i in 1:nrow(table1)){
    if(table6[i,9]!="NA"){
      t6r8=t6r8+1
    }
  }
  r1<-c(t1r1,t2r1,t3r1,t4r1,t5r1,t6r1)
  r2<-c(t1r2,t2r2,t3r2,t4r2,t5r2,t6r2)
  r3<-c(t1r3,t2r3,t3r3,t4r3,t5r3,t6r3)
  r4<-c(t1r4,t2r4,t3r4,t4r4,t5r4,t6r4)
  r5<-c(t1r5,t2r5,t3r5,t4r5,t5r5,t6r5)
  r6<-c(t1r6,t2r6,t3r6,t4r6,t5r6,t6r6)
  r7<-c(t1r7,t2r7,t3r7,t4r7,t5r7,t6r7)
  r8<-c(t1r8,t2r8,t3r8,t4r8,t5r8,t6r8)
  counts<-cbind(r1,r2,r3,r4,r5,r6,r7,r8)
  barplot(counts, xlab="Taxonomic Rank & Table Names", ylab = "Number of Mapped", main="Taxonomy Classification for Ranks",beside=TRUE,col=c("grey0","grey17","grey34","grey51","grey68","grey85"))
  legend("topright",legend=c("idtaxsilva","idtaxpr2","LCAPr2","LCASilva","BayesSilva","BayesPr2"),bty='n',fill=c("grey0","grey17","grey34","grey51","grey68","grey85"))
}

#6 way comparisons: Six-way NA comparison function is created here to count the number of NA for each sequence in each table. A csv file is then created to record. 
sixway<-function(table1,table2,table3,table4,table5,table6){
  countNA<-c()
  for(i in 1:nrow(table1)){
    a=0
    b=0
    c=0
    d=0
    e=0
    f=0
    for(j in 2:ncol(table1)){
      if (table1[i,j]=="NA"){
        a=a+1
      }
      if(table2[i,j]=="NA"){
        b=b+1
      }
      if(table3[i,j]=="NA"){
        c=c+1
      }
      if(table4[i,j]=="NA"){
        d=d+1
      }
      if(table5[i,j]=="NA"){
        e=e+1
      }
      if(table6[i,j]=="NA"){
        f=f+1
      }
    } 
    countNA<-append(countNA,a)
    countNA<-append(countNA,b)
    countNA<-append(countNA,c)
    countNA<-append(countNA,d)
    countNA<-append(countNA,e)
    countNA<-append(countNA,f)
  }
  comparison<-data.frame(matrix(unlist(countNA), nrow=24126, byrow=TRUE))
  colnames(comparison)<-c("Idtaxsilva","Idtaxpr2","LCAPr2","LCASilva","BayesianSilva","BayesianPr2")
  write.csv(comparison,"/Users/cl/Desktop/SixWayComparison.csv")
  View(comparison)
}


AllAgree<-function(table1,table2,table3,table4,table5,table6){
  count=0
  for(i in 1:nrow(table1){
    if (identical(table1[i,],table2[i,])&identical(table1[i,],table3[i,])&identical(table1[i,],table4[i,])&identical(table1[i,],table5[i,])&identical(table1[i,],table6[i,])){
      count=count+1
    }
  }
  print(count)
}

#Three Way Tables Comparison: A three way comparison function is created here to record number of sequences that agree for all taxonomy tanks, for only one, two, and all the way up to seven agreements. This function is used twice to compare all three tables mapped with pr2 and all three tables mapped with silva.

threeway<-function(table1,table2,table3){
  All=0
  One=0
  Two=0
  Three=0
  Four=0
  Five=0
  Six=0
  Seven=0
  for(i in 1:nrow(table1)){
    if(identical(table1[i,],table2[i,])&identical(table2[i,],table3[i,])){
      All=All+1
    }
  }
  for(i in 1:nrow(table1)){
    No=0
    for(j in 2:ncol(table1)){
      if(identical(table1[i,j],table2[i,j])&identical(table2[i,j],table3[i,j])){
        No=No+1
      }
    }
    if(No==1){
      One=One+1
    }
    if(No==2){
      Two=Two+1
    }
    if(No==3){
      Three=Three+1
    }
    if(No==4){
      Four=Four+1
    }
    if(No==5){
      Five=Five+1
    }
    if(No==6){
      Six=Six+1
    }
    if(No==7){
      Seven=Seven+1
    }
  }
  print(All)
  print(One)
  print(Two)
  print(Three)
  print(Four)
  print(Five)
  print(Six)
  print(Seven)
  
  bars1<-c(All, One, Two, Three, Four, Five, Six, Seven)
  lbls1<-c("All Agree","One","Two","Three","Four","Five","Six","Seven","None")
  barplot(bars1, xlab="Number of Agreements",names.arg=lbls1, ylab = "Number of Sequences", col=rainbow(length(bars1)), main="Agreements Among Silva")

}
#Same as the threeway, but with six tables
SixWay<-function(table1,table2,table3,table4,table5,table6){
  None=0
  All=0
  One=0
  Two=0
  Three=0
  Four=0
  Five=0
  Six=0
  Seven=0
  for(i in 1:nrow(table1)){
    if(identical(table1[i,],table2[i,])&identical(table2[i,],table3[i,])&identical(table3[i,],table4[i,])&identical(table4[i,],table5[i,])&identical(table5[i,],table6[i,])){
      All=All+1
    }
  }
  for(i in 1:nrow(table1)){
    No=0
    for(j in 2:ncol(table1)){
      if(identical(table1[i,j],table2[i,j])&identical(table2[i,j],table3[i,j])&identical(table3[i,j],table4[i,j])&identical(table4[i,j],table5[i,j])&identical(table5[i,j],table6[i,j])){
        No=No+1
      }
    }
    if(No==0){
      None=None+1
    }
    if(No==1){
      One=One+1
    }
    if(No==2){
      Two=Two+1
    }
    if(No==3){
      Three=Three+1
    }
    if(No==4){
      Four=Four+1
    }
    if(No==5){
      Five=Five+1
    }
    if(No==6){
      Six=Six+1
    }
    if(No==7){
      Seven=Seven+1
    }
  }
  print(None)
  print(All)
  print(One)
  print(Two)
  print(Three)
  print(Four)
  print(Five)
  print(Six)
  print(Seven)
  
  bars1<-c(All,One,Two,Three,Four,Five,Six,Seven)
  lbls1<-c("All Agree","One","Two","Three","Four","Five","Six","Seven","None")
  barplot(bars1, xlab="Number of Agreements",names.arg=lbls1, ylab = "Number of Sequences", col=rainbow(length(bars1)), main="Agreements Among Six Tables")
}

#confidence boxplot: plot the confidence distribution of the original data
confidence<-c()
for(i in 1:nrow(data1){
  confidence<-append(confidence, data1[[i]]$confidence)
}
confidence1<-c()
for(i in 1:nrow(data2){
  confidence1<-append(confidence1, data2[[i]]$confidence)
}
confidence2<-c()
for(i in 1:nrow(data3){
  for(j in 1:ncol(data3){
    confidence2<-append(confidence2, data3[["boot"]][i,j])
  }
}
confidence3<-c()
for(i in 1:nrow(data4){
  for(j in 1:ncol(data4){
    confidence3<-append(confidence3, data4[["boot"]][i,j])
  }
}

boxplot(confidence,confidence1,confidence3,confidence2,
        main = "Confidence Level Comparison",
        names = c("idtax_silva", "idtax_pr2", "bayes_silva", "bayes_pr2"),
        col = c("orange","red")
)
