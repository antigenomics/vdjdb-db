FROM ubuntu:18.04

# Fix certificate issues
RUN apt-get update && \
    apt-get install -y ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# `openjdk-8`
RUN apt-get update
RUN apt-get install -y --no-install-recommends software-properties-common
# RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EB9B1D8886F44E2A
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

# gradle version
# Gradle 5.1.1!

RUN apt-get install -y unzip

RUN mkdir /opt/gradle
RUN wget https://services.gradle.org/distributions/gradle-5.1.1-all.zip
RUN unzip -d /opt/gradle gradle-5.1.1-all.zip
ENV PATH="${PATH}:/opt/gradle/gradle-5.1.1/bin"
RUN rm gradle-5.1.1-all.zip

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
RUN Rscript -e 'install.packages(c("rmarkdown", "ggplot2", "knitr", "ggpubr", "RColorBrewer", "data.table", "forcats", "ggh4x", "ggalluvial", "ggrepel", "tidyverse", "dplyr", "httr", "xml2", "stringr", "gridExtra"), repos = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'
RUN Rscript -e 'install.packages("reshape2", repose = c("http://cran.us.r-project.org", "https://cloud.r-project.org/"))'

RUN touch .RProfile
RUN echo 'library(rmarkdown)' >> .RProfile
RUN echo 'library(ggplot2)' >> .RProfile
RUN echo 'library(knitr)' >> .RProfile
RUN echo 'library(ggpubr)' >> .RProfile
RUN echo 'library(RColorBrewer)' >> .RProfile
RUN echo 'library(data.table)' >> .RProfile
RUN echo 'library(forcats)' >> .RProfile
RUN echo 'library(ggh4x)' >> .RProfile
RUN echo 'library(ggalluvial)' >> .RProfile
RUN echo 'library(ggrepel)' >> .RProfile
RUN echo 'library(tidyverse)' >> .RProfile
RUN echo 'library(dplyr)' >> .RProfile
RUN echo 'library(httr)' >> .RProfile #
RUN echo 'library(xml2)' >> .RProfile
RUN echo 'library(stringr)' >> .RProfile
RUN echo 'library(gridExtra)' >> .RProfile
RUN echo 'library(reshape2)' >> .RProfile

COPY vdjdb-motifs .

RUN touch docker.sh
RUN echo '# /bin/sh' >> docker.sh
RUN echo '[[ -s "$HOME/.sdkman/bin/sdkman-init.sh" ]] && source "$HOME/.sdkman/bin/sdkman-init.sh"' >> docker.sh
RUN echo 'git clone https://github.com/antigenomics/vdjdb-db' >> docker.sh
RUN echo 'cd vdjdb-db/' >> docker.sh
RUN echo 'bash release.sh' >> docker.sh
RUN echo 'mkdir -p /root/output' >> docker.sh
RUN echo 'cp -r * /root/output/' >> docker.sh

CMD [ "bash", "docker.sh" ]

