.PHONY : all clean build run-macos run-linux

bin_name = font-demo
build_path = build

all: | build

clean:
	@rm -rf $(build_path)

build:
	@mkdir -p $(build_path)
	@cd $(build_path) && cmake .. && make -j4;

run-macos: build
	@cd $(build_path) && ./$(bin_name) "/Library/Fonts/Times New Roman.ttf";

run-linux: build
	@cd $(build_path) && ./$(bin_name) "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf";
