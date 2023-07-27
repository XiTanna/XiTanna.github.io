---
layout: archive
title: "Blog Posts"
permalink: /blogposts/
author_profile: true
---

{% include base_path %}
{% for post in site.blogposts %}
  {% include archive-single.html %}
{% endfor %}