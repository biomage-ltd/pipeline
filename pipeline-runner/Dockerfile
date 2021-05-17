# Create builder step
# pull official base image and use it as builder
FROM r-base:4.0.2 AS builder
WORKDIR /src

# install required debian packages for R installation
ENV DEBIAN_FRONTEND noninteractive
COPY requirements_debian.txt .
RUN apt-get update && apt-get -y upgrade
RUN cat requirements_debian.txt | xargs apt-get -y install --no-install-recommends
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# add GitHub PAT if required for GitHub installations.
ARG GITHUB_PAT
ENV GITHUB_PAT ${GITHUB_PAT}

# install pak for more streamlined installation of packages
RUN R -e 'utils::install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")'

# install required packages
RUN R -e 'pak::pkg_install("remoter")'
RUN R -e 'pak::pkg_install("devtools")'
RUN R -e 'pak::pkg_install("BiocGenerics")'
RUN R -e 'pak::pkg_install("RJSONIO")'
RUN R -e 'devtools::install_version("spatstat", version = "1.64-1")'
RUN R -e 'pak::pkg_install("Seurat")'
RUN R -e 'pak::pkg_install("gprofiler2")'
RUN R -e 'pak::pkg_install("sccore")'
RUN R -e 'pak::pkg_install("GO.db")'

# install handy development stuff
RUN R -e 'pak::pkg_install("zeallot")'
RUN R -e 'pak::pkg_install("paws")'
RUN R -e 'pak::pkg_install("import")'
RUN R -e 'pak::pkg_install("ids")'

# default BLAS/LAPACK issue: Seurat::FindIntegrationAnchors finds ~30 vs ~13,000 anchors
# eventually move this up with other system installs
RUN apt-get update && apt-get -y install libatlas3-base
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# add jq and aws cli so we can get identities using AWS IRSA roles
#
# see at https://docs.aws.amazon.com/eks/latest/userguide/iam-roles-for-service-accounts.html
# for what IRSA roles are
#
# see at https://github.com/paws-r/paws/blob/main/docs/credentials.md for how PAWS
# handles this
#
RUN wget -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 && chmod +x ./jq && cp jq /usr/bin
RUN pip install awscli
COPY aws_config /root/.aws/config

# start app
COPY src ./
ENTRYPOINT ["Rscript", "init.r"]


# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

RUN R -e 'pak::pkg_install("languageserver")'
RUN R -e 'pak::pkg_install("lintr")'
RUN R -e 'pak::pkg_install("styler")'

# start app
COPY src ./
ENTRYPOINT ["Rscript", "init.r"]
