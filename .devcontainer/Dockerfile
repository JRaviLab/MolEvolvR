FROM ghcr.io/rocker-org/devcontainer/r-ver:4

WORKDIR /workspaces/MolEvolvR

# Install MolEvolvR Package Dependencies
## Copy DESCRIPTION to leverage build cache
COPY DESCRIPTION .
RUN R -q -e "pak::local_install_deps('.', dependencies = TRUE)"

## Then copy everything else
COPY . .

CMD ["bash"]
