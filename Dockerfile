FROM r-base:4.5.0

RUN apt-get update
RUN apt-get install -y libgmp-dev libfontconfig1-dev libfreetype6-dev

RUN R -e "install.packages('remotes')"

RUN R -e "remotes::install_github('qBioTurin/Insite')"