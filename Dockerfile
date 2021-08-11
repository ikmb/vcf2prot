FROM docker/whalesay:latest
LABEL Name=ppg Version=0.1.3
RUN apt-get -y update && apt-get install -y fortunes
CMD ["sh", "-c", "/usr/games/fortune -a | cowsay"]
