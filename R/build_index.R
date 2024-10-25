


refs <- dir(system.file(package="Rbowtie2", "extdata", "bt2","refs"),full=TRUE)
(cmdout<-bowtie2_build(references=refs, 
                       bt2Index=file.path(td, "lambda_virus"), "--threads 4 --quiet",
                       overwrite=TRUE))