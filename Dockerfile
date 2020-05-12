FROM rocker/verse:latest

RUN R -e 'install.packages(c("data.table", "raster"))'
RUN R -e 'install.packages(c("testtthat"))'
RUN R -e 'install.packages(c("devtools"))'

COPY . /nosoi

CMD R -e 'devtools::check("nosoi", force_suggests = FALSE, run_dont_test = TRUE)'
