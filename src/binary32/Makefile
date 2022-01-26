FUNCTIONS := acos acosh asin asinh atan atanh

all: $(addprefix check-,$(FUNCTIONS))

check-%:
	$(MAKE) -C $* check

clean:
	for u in $(FUNCTIONS); do $(MAKE) -C $$u clean; done
