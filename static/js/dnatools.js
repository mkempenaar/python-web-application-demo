/**
 * Javascript functions and events
 */

$(function () {
    /**
     * Perform a request for translated sequences and display as a DataTable
     */
    $("#dna-translate-submit-update").click(function (event) {
        // Make sure we don't leave this page
        event.preventDefault();

        /* Create a new FormData object for POSTing to the backend. This will be filled with the
        user-entered data from the form. */
        let form_data = new FormData();
        form_data.append('seq', $("#seq").val());
        
        // Get the HTML elements that are checkboxes
        let checkboxes = document.querySelectorAll('input[type=checkbox]');
        
        /* Add the selected frames to the form-data (using filter and map). Note that the 'checkboxes' object
        is a 'NodeList' which needs to be converted to an Array for iteration. */
        Array.from(checkboxes).filter(checkbox => {
            return checkbox.checked;
        }).map(checkbox => {
            if (checkbox.name == 'fframe')
                form_data.append('fframe', checkbox.value);
            else
                form_data.append('rframe', checkbox.value);
        });

        // TODO: check if both a sequence is entered and frames have been selected. Inform the user if data is missing.

        // Get the translated sequences
        fetch('/dna-translate-data', {
            method: 'POST',
            body: form_data
        })
            .then(function (response) {
                return response.json();
            })
            .then(function (table_data) {
                /* Place the received data into a DataTable (table placeholder is in the HTML).
                A number of default options have been disabled, refer to the DataTables manual for details. */
                $('#translated-table').DataTable({
                    destroy: true,
                    data: table_data[0],
                    columns: [
                        {title: "Frame"},
                        {title: "Sequence"}
                    ],
                    paging: false,
                    searching: false,
                    ordering: false
                })
            });
    });

    /* Perform a request for a GC-percentage Bokeh plot for sequences in an uploaded file */
    $("#gc-submit").click(function (event) {
        // Make sure the form doesn't actually submit
        event.preventDefault();

        // Get the selected file
        let form_data = new FormData();
        let file_field = document.querySelector('input[type="file"]');
        form_data.append('file', file_field.files[0]);

        // TODO: check if a file is selected, ignore if not

        // Get the plot data given the uploaded file
        fetch('/gc-percentage', {
            method: 'POST',
            body: form_data
        })
            .then(function (response) {
                return response.json();
            })
            .then(function (plot) {
                // Empty the div (otherwise plots are appended)
                $("#GCplot").html("");
                // Embed the plot in the given div
                Bokeh.embed.embed_item(plot, "GCplot");
            });
    });
});