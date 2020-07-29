FROM mcr.microsoft.com/dotnet/nightly/sdk:5.0

WORKDIR /fslinalg

# Copy everything and build
COPY . ./
RUN dotnet publish -c Release -o out

# Run tests
RUN dotnet run --project FsLinAlg.Test

ENTRYPOINT ["dotnet", "fsi", "--langversion:preview"]