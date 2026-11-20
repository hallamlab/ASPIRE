document.addEventListener("DOMContentLoaded", () => {
    for (const block of document.querySelectorAll(".highlight-bash, .highlight-shell")) {
        const code = block.querySelector("pre");
        if (!code) continue;

        const button = document.createElement("button");
        button.type = "button";
        button.className = "copy-commands";
        button.textContent = "Copy";
        button.setAttribute("aria-label", "Copy commands");
        button.addEventListener("click", async () => {
            try {
                await navigator.clipboard.writeText(code.textContent.trimEnd());
                button.textContent = "Copied";
            } catch {
                const selection = window.getSelection();
                const range = document.createRange();
                range.selectNodeContents(code);
                selection.removeAllRanges();
                selection.addRange(range);
                button.textContent = "Selected — copy manually";
            }
        });
        button.addEventListener("blur", () => { button.textContent = "Copy"; });
        block.prepend(button);
    }
});
