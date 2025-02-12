const log = function (title, content) {
  let text =
    title +
    " ".repeat(Math.max(log.title_length - title.length, 5)) +
    ":" +
    " ".repeat(log.content_space) +
    content;

  console.log(text);
};

log.title_length = 15;
log.content_space = 5;

export default log;
