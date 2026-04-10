FROM julia:1.11-bookworm

ENV JULIA_DEPOT_PATH=/opt/julia-depot

ENV DEBIAN_FRONTEND=noninteractive

RUN rm /var/lib/dpkg/info/libc-bin.* \
    && apt-get clean \
    && apt-get update \
    && apt-get install -y \
        libcurl4-openssl-dev \
        bzip2 \
        libbz2-dev \
        liblzma-dev \
        git \
        cmake \
        wget \
        procps \
        clang \
        libblas-dev \
        liblapack-dev \
        libopenblas-dev \
        tabix \
        unzip \
        default-jdk \
        zlib1g-dev

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 \
    && tar -xvjf samtools-1.22.1.tar.bz2 \
    && cd samtools-1.22.1 \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Install BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -xvjf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip \
    && unzip plink2_linux_x86_64_20250129.zip \
    && mv plink2 /usr/local/bin \
    && chmod +x /usr/local/bin/plink2

# Install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
    && unzip plink_linux_x86_64_20241022.zip \
    && mv plink /usr/local/bin \
    && chmod +x /usr/local/bin/plink

# Install conda
RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3.sh -b -p /opt/miniforge3 \
    && /opt/miniforge3/bin/conda create -y -n liftover_env python=3 \
    && /opt/miniforge3/bin/conda run -n liftover_env pip install pyliftover pandas \
    && /opt/miniforge3/bin/conda clean -a -y \
    && chmod -R 777 /opt/miniforge3/

COPY . /opt/decodeme

RUN cd /opt/decodeme \
    && julia --project=/opt/decodeme --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' \
    && julia --project=/opt/decodeme --startup-file=no --threads=auto -e 'using PackageCompiler; create_sysimage(["CSV", "DataFrames"], sysimage_path="/opt/decodeme/sysimage.so", include_transitive_dependencies=false, cpu_target=PackageCompiler.default_app_cpu_target())'
