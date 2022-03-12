FROM ubuntu:18.04

# Fix certificate issues
RUN apt-get update && \
    apt-get install -y ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# `openjdk-8`
RUN apt-get update
RUN apt-get install -y --no-install-recommends software-properties-common
RUN add-apt-repository -y ppa:openjdk-r/ppa
RUN apt-get update

RUN apt-get update
RUN add-apt-repository -y universe
RUN apt-get update
RUN apt-get install -y python python-pip

RUN pip install 'pandas==0.24.2' 'numpy==1.16.6'
RUN pip install 'biopython==1.76'

RUN apt-get install -y openjdk-8-jre openjdk-8-jdk openjdk-8-jdk-headless openjdk-8-jre-headless
RUN update-alternatives --config java
RUN update-alternatives --config javac

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

RUN apt-get install -y wget
RUN apt-get install -y unzip
RUN apt-get install -y git
RUN apt-get install -y curl
RUN apt-get install -y zip
RUN apt-get install -y pandoc

SHELL ["/bin/bash", "-c"] 

RUN curl -s "https://get.sdkman.io" | bash
RUN source "$HOME/.sdkman/bin/sdkman-init.sh" && \
    sdk install groovy 3.0.9

# needed for R 
ENV DEBIAN_FRONTEND noninteractive

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
         r-base \
         r-base-dev \
         r-recommended

# for R deps
RUN apt-get install -y libnlopt-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libxml2 libxml2-dev

# 'knitr', 'htmltools', 'jquerylib', 'stringr' are not available for package 'rmarkdown'
RUN Rscript -e 'install.packages(c("rmarkdown", "ggplot2", "knitr", "ggpubr", "RColorBrewer", "data.table", "forcats", "ggh4x", "ggalluvial", "ggrepel", "tidyverse", "dplyr", "httr", "xml2", "stringr", "gridExtra", "circlize", "maps", "scatterpie"), repos = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'
RUN Rscript -e 'install.packages("reshape2", repos = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'
RUN Rscript -e 'install.packages(c("stringdist", "ggseqlogo", "igraph"), repos = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'
RUN Rscript -e 'install.packages(c("reshape2", "FField", "reshape", "gplots", "gridExtra", "circlize", "ggplot2", "grid", "VennDiagram", "ape", "MASS", "plotrix", "RColorBrewer", "scales"), repos = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'

RUN apt-get install -y texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
RUN apt-get install -y build-essential procps curl file git

RUN wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
RUN mkdir -p /software/bin/
RUN unzip vdjtools-1.2.1.zip
RUN cp -r vdjtools-1.2.1/* /software/bin/
RUN chmod +x /software/bin/vdjtools

ENV PATH="/usr/lib/jvm/java-8-openjdk-amd64/bin:${PATH}"
ENV PATH="/software/bin:${PATH}"

COPY vdjdb-motifs vdjdb-motifs

RUN touch docker.sh
RUN echo '# /bin/sh' >> docker.sh
RUN echo '[[ -s "$HOME/.sdkman/bin/sdkman-init.sh" ]] && source "$HOME/.sdkman/bin/sdkman-init.sh"' >> docker.sh
RUN echo 'git clone https://github.com/antigenomics/vdjdb-db' >> docker.sh
RUN echo 'cd vdjdb-db/' >> docker.sh
RUN echo 'bash release.sh' >> docker.sh
RUN echo 'mkdir -p /root/output' >> docker.sh
RUN echo 'cp -r database/*zip /root/output/' >> docker.sh

CMD [ "bash", "docker.sh" ]