From noaaepic/ubuntu20.04-intel-landda:develop

CMD ["/bin/bash"]

ENV HOME=/opt
WORKDIR $HOME

#remove org land-offline_workflow
RUN rm -rf $HOME/land-offline_workflow
COPY . $HOME/land-offline_workflow

# set env vars
ENV FIXlandda=$HOME/land-DA_workflow/fix
ENV EPICHOME=/opt
ENV JEDI_INSTALL=${EPICHOME}/jedi_skylabv7.0
ENV TEST_BASEDIR=${EPICHOME}/test_base/restarts/vector"}

#build & unit testing
WORKDIR $HOME/land-offline_workflow
RUN source /opt/spack-stack/.bashenv; mkdir build; cd build; pwd; ecbuild ..; make -j2; ctest -V --stop-on-failure
