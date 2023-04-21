PREFIX = /usr/share/asymptote
TARGET = jomatra.asy
MODULES = geo-modules

install:
	install -m 644 $(TARGET) $(PREFIX)
	install -d $(MODULES) $(PREFIX)

uninstall:
	rm -rf $(PREFIX)/$(MODULES) $(PREFIX)/$(TARGET)

test:
	find tests -type f -name "*.asy" -exec asy -dir . -o ./build {} \;

clean:
	find . -type f -name "*.eps" -exec rm {} \;
	find . -type f -name "*.pdf" -exec rm {} \;
	rm -r build
