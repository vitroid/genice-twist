.DELETE_ON_ERROR:
GENICE=genice2
BASE=genice2_twist
PACKAGE=genice2-twist

all: README.md

%: temp_% replacer.py $(wildcard $(BASE)/formats/*.py) $(BASE)/__init__.py
	pip install genice2_dev
	python replacer.py < $< > $@

test: CRN1.btwc.test CRN1.CM.twist.png CRN1.DB.twist.png CRN1.SB.twist.png CRN1.twist.png CRN1.twist.svg\
      CRN1.twist.yap # CRN1.CM.twist.svg CRN1.DB.twist.svg CRN1.SB.twist.svg

%.test:
	make $*
	diff $* ref/$*
%.btwc: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist > $@
%.CM.twist.png: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[png:CM] > $@
%.CM.twist.svg: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[svg:CM] > $@
%.DB.twist.png: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[png:DB] > $@
%.DB.twist.svg: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[svg:DB] > $@
%.SB.twist.png: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[png:SB] > $@
%.SB.twist.svg: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[svg:SB] > $@
%.twist.yap: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[yaplot] > $@
%.twist.png: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[png] > $@
%.twist.svg: $(BASE)/formats/twist.py
	$(GENICE) $* -f twist[svg] > $@
sample.png: CRN1.twist.svg
	inkscape $< -o $@

test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PACKAGE)



install:
	python ./setup.py install
uninstall:
	-pip uninstall -y $(PACKAGE)
build: README.md $(wildcard $(BASE)/formats/*.py)
	python ./setup.py sdist bdist_wheel


deploy: build
	twine upload --repository pypi dist/*
check:
	./setup.py check
clean:
	-rm $(ALL) *~ */*~ *.btwc *.svg *.yap CRN1*.png ja
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
