library(data.table)
library(ggplot2)
library(dplyr)
#library(ggpubr)
#library(ggsci)
library(stringr)
library(tidyr)

control_medication_measures = c("statin_LDL", "statin_HDL", "statin_TC", "metformin_HbA1c", "betablocker_HR", 
                               "antihpt_all_SBP") # traits for medication-free individuals

filter = "l"
n_meas = "1"

for (pair in control_medication_measures){

    for (out_pheno in c("pheno_diff", "pheno_log_diff", "pheno_prop", "pheno_diff_adj")){

        if (out_pheno == "pheno_diff"){
            pheno_name = "post-base"
        } else if (out_pheno == "pheno_log_diff"){
            pheno_name = "log(post)-log(base)"
        } else if (out_pheno == "pheno_prop"){
            pheno_name = "post/base"
        } else if (out_pheno == "pheno_diff_adj"){
            pheno_name = "post-base, base adj."
        }

        X_all_sub <- fread(paste0("output/GWAS/control_", pair, "_", filter, "_", n_meas, "_GWASformatted_", out_pheno, ".tsv"))
        X_all_sub$LOG10P = -log10(X_all_sub$p)
        X_all_sub = X_all_sub[X_all_sub$LOG10P > 1]
        
        if (str_detect(pair, "statin")){
            drug = "statin"
            phenotype = str_split(pair, "_")[[1]][2]
            my_colors <- c("#01255e", "#505c70")
        } else if (str_detect(pair, "metformin")){
            drug = "metformin"
            phenotype = "HbA1c"
            my_colors <- c("#7d1a15", "#7d6447")
        } else if (str_detect(pair, "antihpt")){
            phenotype = "SBP"
            drug = str_split(pair, "_")[[1]][2]
            if (drug == "thiazide"){
                drug = "Thiazide d."
            } else if (drug == "all"){
                drug = "All antih."
            }
            my_colors <- c("#025208", "#6e803d")
        } else if (str_detect(pair, "betablocker")){
            drug = "beta bl."
            my_colors <- c("#63043f", "#694a5d")
            if (str_detect(pair, "SBP")){
                phenotype = "SBP"
            } else if (str_detect(pair, "HR")){
                phenotype = "HR"
            } 
        } 
        
        N = X_all_sub[1, "N"]
        
        title <- paste0(phenotype, "-longitudinal control, ", pheno_name, "; N = ", N)
        
        X_all_sub[,max_bp:=max(GENPOS),by="CHROM"]
        toAdd <- unique(X_all_sub[,list(CHROM,max_bp)])
        toAdd[,toAdd:=lag(max_bp)]
        toAdd[is.na(toAdd),toAdd:=0]
        toAdd[,cumToAdd:=cumsum(as.numeric(toAdd))]
        
        X_all_sub <- merge(X_all_sub,toAdd[,list(CHROM,cumToAdd)],by="CHROM",all.x=T)
        X_all_sub[,bpAdj:=GENPOS+cumToAdd]
        
        axis_set <- X_all_sub %>% 
            group_by(CHROM) %>% 
            summarize(center = mean(bpAdj))
        
        X_all_sub[,CHR2:=" "]
        X_all_sub[(CHROM)%%2 == 0,CHR2:="  "]
        
        axis_set_mod <- as.data.table(axis_set)
        axis_set_mod[,CHR_write:=as.character(CHROM)]
        axis_set_mod <- axis_set_mod[!(CHROM %in% c(11,13,15,16,18,19,21,22)),]
        axis_set_mod[CHR_write == "23",CHR_write:="X"] 
        
        
        X_all_sub_ordered <- copy(X_all_sub)
        
        X_all_sub_ordered[LOG10P>249,LOG10P:=249]
        
        # add names so bar for 'c' gets fill, too
        names(my_colors) <- c(" ","  ")
        p_manhattan <- ggplot(data=X_all_sub_ordered, aes(x=bpAdj,y=LOG10P,col=CHR2)) + geom_point(size=1) + 
            theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2,'cm')) +
            #theme(legend.title = element_blank(), legend.text = element_text(size=13), legend.key.size = unit(0.5,'cm'), legend.spacing.x = unit(0.2,'cm')) +
            theme(title = element_text(color="gray20", size=13)) +
            theme(axis.title = element_text(color="gray20", size=13)) +
            theme(axis.text = element_text(color="gray40",size=13)) +
            #theme(legend.position = "top") +
            theme(strip.text=element_text(size = 13, color = "gray10")) +
            theme(panel.background = element_blank(), rect = element_rect(fill = "transparent")) +
            xlab("Chromosome") + ylab("-log10(p)") + geom_hline(col="black",yintercept=-log10(5e-8),lty=3) +
            scale_color_manual(values = my_colors) + #,breaks=c("Novel discovery","","") #,"col_odd","col_even"
            scale_x_continuous(label = axis_set_mod$CHR_write, breaks = axis_set_mod$center) + 
            ylim(c(0,max(X_all_sub_ordered$LOG10P)+5)) + guides(color = "none") + 
            ggtitle(title)
        
        png(paste0("output/Manhattan/control_", pair, "_", filter, "_", n_meas, "_GWAS_", out_pheno, "_manhattan.png"), width = 8, height = 5, res = 300, units = "in")
        print(p_manhattan)
        dev.off()
    }
}

