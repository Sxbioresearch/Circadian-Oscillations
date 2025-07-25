###Difference analysis
##Absolute quantification
ZicoSeq1 <- ZicoSeq(meta.dat = meta, feature.dat = as.matrix(input), 
                    grp.name = 'group',feature.dat.type = "other",
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct=0.01,winsor.end = 'both',
                    # Posterior sampling 
                    is.post.sample = F,
                    # Use the square-root transformation
                    stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata =  meta$Env , 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
##relative abundance
ZicoSeq1 <- ZicoSeq(meta.dat = meta_r, feature.dat = as.matrix(input_r), 
                    grp.name = 'group',feature.dat.type = "proportion",
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct=0.01,winsor.end = 'top',
                    # Posterior sampling 
                    is.post.sample = TRUE,
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5),stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata =  meta_r$Env , 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
##foldchange
foldChange<-function(inData,classLabel)
{
  #Calculating all probes' FC value
  sampleIdsCase<-which(classLabel=="timepoint1");
  sampleIdsControl<-which(classLabel=="timepoint2");
  probeFC<-rep(0,nrow(inData))
  for(i in 1:nrow(inData))
  {
    probeFC[i]<-mean(as.numeric(as.vector(inData[i,sampleIdsCase])))/mean(as.numeric(as.vector(inData[i,sampleIdsControl])));
  }
  probeFC<-log(probeFC,base=2);
  result<-probeFC; 
  return(result)
}


