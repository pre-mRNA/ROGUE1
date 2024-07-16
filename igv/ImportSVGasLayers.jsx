// Batch SVG Import - Adobe Illustrator Script
// Purpose: This script imports SVG files from a selected folder into a new Illustrator document as individual layers,
// scales them down to 3% of their original size, crops to visible bounds, and aligns them to the top left corner of an A4 artboard.
// Version: 1.0.4 - 2024-07-16

function selectImportFolder() {
    return Folder.selectDialog('Choose the folder containing the SVG files:');
}

function batchImportSVG(selectedFolder) {
    var doc;
    if (selectedFolder) {
        // Define A4 dimensions in points (1 point = 1/72 inch)
        var a4Width = 595.28; // 210mm in points
        var a4Height = 841.89; // 297mm in points
        doc = app.documents.add(DocumentColorSpace.RGB, a4Width, a4Height);

        var files = selectedFolder.getFiles(function(f) {
            return f instanceof File && f.name.toLowerCase().search(".svg") > -1;
        });

        for (var i = 0; i < files.length; i++) {
            var file = files[i];
            var layer = doc.layers.add({ name: file.name.split('.')[0] });
            var svgItem = layer.groupItems.createFromFile(file);

            svgItem.resize(3, 3, true, true, true, true, 3, Transformation.DOCUMENTORIGIN);
            clearClippingMasks(svgItem);
            fitArtworkToBounds(svgItem);

            // Move each SVG to the top left corner of the artboard
            svgItem.position = [0, a4Height - svgItem.height];
        }

        if (doc.layers.length === 0) {
            alert("No SVG files found. Please verify the folder and try again.");
            doc.close();
        }
    } else {
        alert("Folder selection was cancelled. Please run the script again.");
    }
}

function clearClippingMasks(item) {
    for (var i = item.pageItems.length - 1; i >= 0; i--) {
        if (item.pageItems[i].clipping) {
            item.pageItems[i].remove();
        }
    }
}

function fitArtworkToBounds(group) {
    var bounds = group.visibleBounds; // [left, top, right, bottom]
    var width = bounds[2] - bounds[0];
    var height = bounds[1] - bounds[3];
    group.left = bounds[0];
    group.top = bounds[1];
    group.width = width;
    group.height = height;
}

// Initialize the script
batchImportSVG(selectImportFolder());
