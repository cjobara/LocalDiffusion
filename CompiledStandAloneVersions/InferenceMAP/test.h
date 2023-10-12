			{
				native.title("Select tr x y t File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_FILE);
				native.filter("Text Files\t*.trxyt");

				// in case save screen is canceled
				if (native.show() == 0) {

					const int fileCount = native.count();

					strcpy(nativeFilename,native.filename());

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}




			{
				native.title("Select xyt File(s)");
				native.type(Fl_Native_File_Chooser::BROWSE_MULTI_FILE);
				native.filter("xyt Files\t*.xyt");
				native.show();
				strcpy(nativeFilename,native.filename());

				const int fileCount = native.count();

				// load file conditions (to avoid closing file if open screen is cancelled)
				if ((nativeFilename[0]!=(char)NULL) || (nativeFilename[0]!=(char)NULL && iMAP->fileNameStorage[iMAP->fileNumber][0]!=(char)NULL)) {

					list = new char*[fileCount];//(char**)calloc(fileCount,sizeof(char*));
					for (int e = 0; e < fileCount; e++) {
						list[e] = new char[FILENAME_MAX];//(char*)calloc(FILENAME_MAX,sizeof(char));
					}

					// copy files to list
					for (int f = 0; f < fileCount; f++) {
						strcpy(list[f],native.filename(f));
						// printf("%s\n",list[f]);
					}

					// store filename
					strcpy(iMAP->fileNameStorage[iMAP->fileNumber], native.filename(0));

					// instantiate plotObject
					iMAP->fileObjectArray[iMAP->fileNumber] = new File(fileCount, list, fileType);

					// refresh all widgets
					Fl::redraw();
					fileOpened = true;
				}
				break;
			}