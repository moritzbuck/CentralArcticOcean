library(data.table)
library(ggplot2)
library(qualpalr)


raw_reads = fread("cat/read_counts_all_samples.csv")

md = read.csv("metadata.csv")
sample_map = read.csv("sample_map.csv")
library_md = merge(sample_map, md, by.x="sample", by.y="X")
lib_sizes = merge(melt(raw_reads), library_md, by.x = "variable", by.y = "X")[, .(lib_sizes = sum(value)), by = .(sample)]
setkey(lib_sizes, "sample")
data = melt(raw_reads, id.vars="V1", variable.name = "library", value.name = "raw_reads")
colnames(data)[1] = "MAB"


mab_md = fread("all_good_mags.stats")
mab_md$type = "coass"
mab_md$type[grepl("P6404",mab_md$bin_name) |  grepl("Sample",mab_md$bin_name)] = "single"

mab_md[,mab_id := paste(type, assmbler, bin_name, sep="_")]
mab_md[,mab_id := gsub("-","_",mab_id)]

setkey(mab_md, mab_id)

tax = fread("sourmash_lca.combo.tax", sep=",")
tax[,MAGs := gsub("-","_",MAGs)]

mab_md = merge(mab_md, tax, by.x = "mab_id", by.y = "MAGs", all.x = TRUE )

viral_dat = fread("viral_fraction.csv"); setkey(viral_dat, mab_id)
mab_md = merge(mab_md, viral_dat, by = "mab_id", all.x = TRUE )
mab_md[, viral_fract := viral_length/length]

data[, MAB := gsub("-","_",MAB)]

data = merge(data, library_md, by.x = "library", by.y = "X")
data$read_fract = data$raw_reads/lib_sizes[as.character(data$sample)]$lib_sizes

data = merge(data, mab_md, by.x = "MAB", by.y = "mab_id", all.x = TRUE)
bla = data[,.(read_fract=sum(read_fract), max_fract = max(read_fract)), by = .(MAB)]
setkey(bla,MAB)
mab_md[,tot_cov := bla[mab_md]$read_fract]
mab_md[,max_cov := bla[mab_md]$max_fract]

data$set = "MAB"
data[!is.na(viral_fract) & viral_fract > 0.1, Phylum := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Class := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Order := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, Family := NA]
data[!is.na(viral_fract) & viral_fract > 0.1, set := "viral" ]
data[MAB == "*", set := "unassembled"]
data[MAB == "unbinned", set := "unbinned"]

pca = predict(prcomp(-log(fract_mat[,-1]+1)))
row.names(pca) = fract_mat[,1]

mds = metaMDS(t(fract_mat[,-1]))$point

dd = as.data.table(merge(as.data.frame(mab_md), as.data.frame(mds), by.x = "mab_id", by.y = "row.names") )


rare_proc = 0.001
temp = data[, .(read_fract = sum(read_fract), fract = fract[1], station = station[1], type = type.x[1]), by = .(Domain, Phylum, Class, Order, Family,Genus, library, set)]

temp = temp[, tax_str := paste(Domain,Phylum, Class,Order,Family, sep = ";" ) ]
temp[set == 'viral', tax_str := "viral"]
temp = temp[set != 'unbinned' & set != 'unassembled' & tax_str != "NA;NA;NA;NA;NA"]
temp[, read_fract := read_fract/sum(read_fract) ,by = .(station, type, fract)]
rares = temp[, .(read_fract = sum(read_fract)),by = .(tax_str, station, type, fract)][, .(max_fract = max(read_fract)), by = tax_str][max_fract < rare_proc]$tax_str
temp[tax_str %in% rares, tax_str := "rares"]

tax_groups = sapply(levels(factor(paste(temp$Domain,temp$Phylum,temp$Class, sep=";"))), function(x) grep(x, levels(factor(temp$tax_str)), val=TRUE) )
tax_groups = tax_groups[sapply(tax_groups, length) > 0 ]
start = qualpal(length(tax_groups), list(h=c(0,360), s = c(0.7,0.7), l=c(0.6, 0.6)))$HSL[,1]
names(start) = names(tax_groups)
cols = lapply(names(start), function(x) if(length(tax_groups[[x]]) > 1) qualpal(n = length(tax_groups[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.7,0.6) )
greys = qualpal(n=2, list(h=c(200,200), s=c(0,0), l=c(0.4, 0.85)))$hex

cols = c(unlist(cols), greys)

ggplot(temp, aes(x = factor(station), y=read_fract, fill = tax_str))+theme_minimal()+ geom_bar(stat = 'identity')+facet_wrap(type~factor(fract), scales = "free")+guides(fill=guide_legend(ncol=2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text=element_text(size=7), legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=cols)


temp = data[, .(read_fract = sum(read_fract), fract = fract[1], station = station[1], type = type.x[1]), by = .(Domain, Phylum, Class, Order, Family,Genus, library, set)]

temp = temp[, tax_str := paste(Domain,Phylum, Class, sep = ";" ) ]
temp$tax_str[temp$tax_str == "NA;NA;NA"] = "unclassified MAB"
temp$tax_str[temp$set == 'unbinned']  = "unbinned"
temp$tax_str[temp$set == 'unassembled']  = "unassembled"
temp$tax_str[temp$set == 'viral']  = "viral"
rares = temp[, .(read_fract = sum(read_fract)),by = .(tax_str, station, type, fract)][, .(max_fract = max(read_fract)), by = tax_str][max_fract < rare_proc]$tax_str
temp[tax_str %in% rares, tax_str := "rares"]

tax_groups = sapply(levels(factor(paste(temp$Domain,temp$Phylum,temp$Class, sep=";"))), function(x) grep(x, levels(factor(temp$tax_str)), val=TRUE) )
tax_groups = tax_groups[sapply(tax_groups, length) > 0 ]

cols = lapply(names(tax_groups), function(x) if(length(tax_groups[[x]]) > 1) qualpal(n = length(tax_groups[[x]]), list(h = c(start[x], start[x]), s = c(0.6, 0.9), l = c(0.6, 0.85)))$hex else hsv(start[[x]]/360, 0.8,0.8) )
greys = qualpal(n=4, list(h=c(200,200), s=c(0,0), l=c(0.4, 0.85)))$hex

cols = c(unlist(cols), rev(greys))

ggplot(temp, aes(x = factor(station), y=read_fract, fill = tax_str))+ geom_bar(stat = 'identity')+theme_minimal()+facet_wrap(type~factor(fract), scales = "free")+guides(fill=guide_legend(ncol=2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text=element_text(size=7), legend.key.size = unit(0.5,"line"))+guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=cols)
