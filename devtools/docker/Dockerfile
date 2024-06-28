FROM molssi/mamba141

COPY ./lammps_step/data/seamm-lammps.yml /root/environment.yml

RUN mamba env update -f /root/environment.yml

WORKDIR /home
ENTRYPOINT ["lmp_serial"]
