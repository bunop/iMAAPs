clust<-function(x){
	y<-x*0;
	y[1]<-1;
	for(i in 2:length(x)){
		if((x[i]-x[i-1]<2)|((x[i]-x[i-1])/(x[i]+x[i-1]))<0.1){
			y[i]<-y[i-1];
		}
		else
		{
			y[i]<-y[i-1]+1;
		}
	}
	return(y);
}

iMAAPsmry<-function(ad_inp, rawld_inp, num_chr, denoise_cut, method="median")
{
# Required arguments:
# ------------------
# ad_inp is a matrix of the output fitting coeffients
# rawld_inp is a matrix of the output Weighted LD decays
# denoise_cut is a number between 0 and 1 define the noise cut threshold
# num_chr is a number of the chromosome used for Jacknife process
# method is a function used to estimate the average of the admixture signal, it has two choices: median and mean

	cutoff<-1-denoise_cut;
	coef_max<-0;
	coef_min<-1;
	#ad_inp<-read.csv(file=ad_inp,header=T,sep="\t");
	#rawld_inp<-read.csv(file=rawld_inp,header=T,sep="\t");
	if(max(rawld_inp$Distance)>5){
		rawld_inp$Distance<-rawld_inp$Distance/100;
	}
	max_dis<-max(rawld_inp$Distance);
	min_dis<-min(rawld_inp$Distance);
#print(c(min_dis, max_dis))

	data<-matrix(data = NA, nrow = dim(ad_inp)[1], ncol = 2);
#data1<-matrix(data = NA, nrow = dim(ad_inp)[1], ncol = 1);
	p_value<-matrix(data = NA, nrow = dim(ad_inp)[1], ncol = 1);
#z_score<-matrix(data = NA, nrow = dim(ad_inp)[1], ncol = 1);
	gen<-ad_inp[,1];
	GEN_index<-1;

##cut off the noise
	da_raw<-0*ad_inp;
	w<-seq(1,dim(ad_inp)[1],1)*0;
	w[1]<-max_dis;
	for(j in 2:dim(ad_inp)[1]){w[j]<-((1-min_dis)^gen[j]-(1-max_dis)^gen[j])/gen[j];}##begin and end site;
	for(k in 2:(2+num_chr)){
		dd<-cbind(seq(1,dim(ad_inp)[1],1), ad_inp[,k]*w);
		sm<-sum(dd[,2]);
		yy<-dd[sort.list(dd[,2],decreasing=T), ];
		sm_new<-sm;
		for( i in 1:dim(ad_inp)[1]){
			sm_new<-sm_new-yy[i,2];
#		if(sm>0){if(((sm_new/sm)>(1-cutoff)))break;}
			if(sm>0){if(((1-sm_new/sm)>(1-cutoff)))break;}

		}
		for(ll in yy[1:i,1]){
			da_raw[ll,k]<-ad_inp[ll,k];#print(ll);
		}

	}
	if(method=="median"){data[,1]<-apply(da_raw[, 3:(2+num_chr)],1, median);}
	if(method=="mean"){data[,1]<-apply(da_raw[, 3:(2+num_chr)],1, mean);}
	data[,2]<-apply(da_raw[, 3:(2+num_chr)],1, sd);
	for(i in 1:dim(ad_inp)[1]){
		p_value[i,GEN_index]<-(t.test(da_raw[i, 3:(2+num_chr)], mu=0, alternative="greater")$p.value);
	}

	slist<-data[,1]>0;
	x<-gen[slist];
	w1<-data[slist,1];
	sum_spectrum<-cbind(x,w1);
	colnames(sum_spectrum)<-c("gen", "mag");
	w2<-w[slist];
	tot<-sum(w1*w2);
	p_data<-p_value[slist,1];
	sd_data<-data[slist,2];
	if(length(x)>=2){
		y<-clust(x);
		nclus<-max(y);
		m<-vector(length=nclus);
		v<-vector(length=nclus);
		p<-vector(length=nclus);
		mag<-vector(length=nclus);
		mag_sd<-vector(length=nclus);
		prop<-vector(length=nclus);
		for (i in 1:nclus){
			m[i]<-round(sum(x[y==i]*w1[y==i])/sum(w1[y==i]),digits = 6);
			v[i]<-round(sqrt(sum((x[y==i]-m[i])^2*w1[y==i]/sum(w1[y==i]))),digits = 6);
			p[i]<-min(p_data[y==i]);
			mag[i]<-max(w1[y==i]);
			mag_sd[i]<-sd_data[w1==max(w1[y==i])];
			prop[i]<-sum((w1*w2)[y==i])/tot;
		}
		time<-cbind(m, v, p, mag, mag_sd, prop);
	}
	if(length(x)==1) time<-cbind(x,0, min(p_data),max(w1), max(sd_data), 1);
	if(length(x)==0) { time<-cbind(NA,NA,NA,NA,NA, NA);}
	colnames(time)<-c(method, "SD", "pvalue", "mag", "mag_sd", "prop");
	outlist<-list();
	outlist[[ "sum_spectrum" ]]<-sum_spectrum;
	outlist[[ "time" ]]<-time;
	return(outlist);

}
