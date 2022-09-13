FROM ubuntu:focal

RUN apt-get update

# Ensure tzdata installation (a dependency) is non-interactive.
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get install -y software-properties-common wget build-essential cmake \
    git python3-dev xxd libgmp-dev nlohmann-json3-dev libssl-dev
RUN apt-get install -y --no-install-recommends libntl-dev

WORKDIR /home
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
RUN tar --bzip2 -xf boost_1_73_0.tar.bz2 && rm boost_1_73_0.tar.bz2
RUN cd /home/boost_1_73_0 && ./bootstrap.sh --with-python=/usr/bin/python3
RUN cd /home/boost_1_73_0 && ./b2 install

WORKDIR /home
RUN git clone https://github.com/emp-toolkit/emp-tool.git \
  && cd emp-tool \
  && cmake -DCMAKE_BUILD_TYPE=Release . \
  && make \
  && make install

WORKDIR /home
RUN git clone https://github.com/emp-toolkit/emp-ot.git \
  && cd emp-ot \
  && cmake -DCMAKE_BUILD_TYPE=Release . \
  && make \
  && make install

CMD ["bash"]
