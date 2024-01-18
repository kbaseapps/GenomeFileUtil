FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------

RUN pip install --upgrade --extra-index-url=https://pypi.anaconda.org/kbase/simple \
  pip \
  biopython==1.70 \
  releng-client==0.0.1 \
  pytest==7.0.1

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
