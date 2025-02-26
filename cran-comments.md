## Build settings

> # Print R version
> version

```
               _                           
platform       x86_64-apple-darwin20       
arch           x86_64                      
os             darwin20                    
system         x86_64, darwin20            
status                                     
major          4                           
minor          4.0                         
year           2024                        
month          04                          
day            24                          
svn rev        86474                       
language       R                           
version.string R version 4.4.0 (2024-04-24)
nickname       Puppy Cup                   
```

> cat("R version:", R.version.string, "\n\n")

```
R version: R version 4.4.0 (2024-04-24) 
```

> # Print full session information
> print(sessionInfo())

```
R version 4.4.0 (2024-04-24)
Platform: x86_64-apple-darwin20
Running under: macOS Monterey 12.6.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/fr_CH.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] archipelago_0.0.0.9000 revdepcheck_1.0.0.9002 dplyr_1.1.4           
[4] ggplot2_3.5.1          testthat_3.2.1.1      
```
> devtools::check_win_release()
> devtools::check_win_devel()
> devtools::check_win_oldrelease()

```
this notification has been generated automatically.
Your package archipelago_0.0.0.9000.tar.gz has been built (if working) and checked for Windows.
Please check the log files and (if working) the binary package at:
https://win-builder.r-project.org/rPBx6fyHGhOY
The files will be removed after roughly 72 hours.
Installation time in seconds: 3
Check time in seconds: 42
Status: 3 WARNINGs, 5 NOTEs
R version 4.3.3 (2024-02-29 ucrt)

All the best,
(CRAN maintainer of binary packages for Windows)
```
