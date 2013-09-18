CHEASE  = chease_src
PROFILE = profile_src
FLUX    = flux_src
THAWC   = tomuhawc_src
GGJ     = ggj_src
STAGE3  = Stage3

MY_TARGETS = $(CHEASE) $(PROFILE) $(FLUX) $(THAWC) $(GGJ) $(STAGE3)

.PHONY: all $(MY_TARGETS)
all: 	$(MY_TARGETS)

$(CHEASE):
	make -C $(CHEASE)

$(PROFILE):
	make -C $(PROFILE)

$(FLUX):
	make -C $(FLUX)

$(THAWC):
	make -C $(THAWC)

$(GGJ):
	make -C $(GGJ)

$(STAGE3):
	make -C $(STAGE3)

clean:
	make -C $(CHEASE)  clean
	make -C $(PROFILE) clean
	make -C $(FLUX)    clean
	make -C $(THAWC)   clean
	make -C $(GGJ)     clean
	make -C $(STAGE3)  clean

clear:
	make -C $(PROFILE) clear
	make -C $(FLUX)    clear
	make -C $(THAWC)   clear
	make -C $(GGJ)     clear
	make -C $(STAGE3)  clear


